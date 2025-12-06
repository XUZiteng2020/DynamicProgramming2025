%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: Huggett_transition
% 
% Author: Kiyea Jin
% Date: Nov 30, 2025
%
% Description:
% This MATLAB script solves the transition dynamics of the Huggett model,
% when there is a permanent increase in unemployment risk.
%
% Reference: huggett_transition.m by Benjamin Moll
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-sigma))/(1-sigma)
% - Elasticity of intertemporal substitution (sigma): 2
% - Discount rate (rho): 0.05
% - Income: z = [z_u, z_e] = [0.1, 0.2];
% - Lambda: la = [la_u, la_e] = [0.6, 0.5] -> [0.6, 0.8];
% - Discrete grid of asset levels (a): -0.15 to 4
% - Discrete grid of time (time): 0 to 20
% - Borrowing constraint: a>=-0.15
% - Delta = 1000; (Can be arbitrarily large in implicit method)
% - the speed of updating the interest rate:
% xi = 20*(exp(-0.05*(1:p.It)) - exp(-0.05*p.It));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1-1. INITIAL EQUILIBRIUM
    
p = define_parameters_initial();
[ss_initial] = Huggett_GE(p);
    
%% 1-2. TERMINAL EQUILIBRIUM
% Notes: The shock is not mean-reverting;
% it represents a one-time change that permanently shifts the parameter. 
% Hence, we use the updated parameters to solve for transition dynamics.

p = define_parameters_terminal(); 
[ss_terminal] = Huggett_GE(p);

%% 2. INITIALIZE GRID POINTS

a = linspace(p.amin, p.amax, p.I)';
da = (p.amax-p.amin)/(p.I-1);

aa = [a, a]; % I*2 matrix
zz = ones(p.I,1).*p.zz; % I*2 matrix

%% TIME

time = linspace(p.tmin, p.tmax, p.It)';
dt = (p.tmax-p.tmin)/(p.It-1);

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the forward and backward finite-difference operators  
% Df such that Df*V=dVf and Db such that Db*V=dVb

    Df = zeros(p.I, p.I);
    Db = zeros(p.I, p.I);

    for i=1:p.I-1
        Df(i,i) = -1/da; Df(i,i+1) = 1/da;
    end
    Df = sparse(Df);

    for i=2:p.I
        Db(i,i-1) = -1/da; Db(i,i) = 1/da;
    end
    Db = sparse(Db);

% 3-2. Construct A matrix

    A = [-p.lambda(1).*speye(p.I), p.lambda(1).*speye(p.I);
         p.lambda(2)*speye(p.I)  , -p.lambda(2)*speye(p.I)]; % 2I*2I matrix

%% TERMINAL CONDITION FOR V
    
    V_t = zeros(p.I, 2, p.It);
    V_T = ss_terminal.V;
    V_t(:,:,end) = V_T;

%% INITIAL CONDITION FOR g

    gg_t = zeros(p.I, 2, p.It);
    gg_0 = ss_initial.gg;
    gg_t(:,:,1) = gg_0;

%% INTEREST RATE ITERATION

% Guess an initial series of interest rates
    
    r_T = ss_terminal.r;
    r_t_0 = r_T*ones(p.It,1);
    r_t = r_t_0; % It*1 column vector

    S_t = zeros(p.It,1);
    dS_t = zeros(p.It,1);

    xi = 20*(exp(-0.05*(1:p.It)) - exp(-0.05*p.It));

for nr=1:p.Nr

    r_t_nr(:,nr) = r_t;

%% 4. SOLVE HJB BACKWARD GIVEN A SERIES OF INTEREST RATES AND TERMINAL CONDITION

for t=p.It:-1:1

    V = V_t(:,:,t);

    % 4-1. Compute the derivative of the value function 
    dVf = Df*V; % I*2 matrix
    dVb = Db*V;

    % 4-2. Boundary conditions
    % a>=a_min is enforced (borrowing constraint)
    dVb(1, :) = p.mu(p.zz + r_t(t)*aa(1,:));
   
    % a<=a_max is enforced which helps stability of the algorithm
    dVf(end, :) = p.mu(p.zz + r_t(t)*aa(end,:));

    I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)

    % 4-3. Compute the optimal consumption
    cf = p.inv_mu(dVf); % I*2 matrix
    cb = p.inv_mu(dVb);
    
    % 4-4. Compute the optimal savings
    sf = zz + r_t(t).*aa - cf; % I*2 matrix
    sb = zz + r_t(t).*aa - cb;

    % 4-5. Upwind scheme
    If = sf>0;  % I*2 matrix
    Ib = sb<0;
    dV0 = p.mu(zz + r_t(t).*aa);

    dV = dVf.*If + dVb.*Ib + dV0.*(1-If-Ib);  % I*2 matrix

    c = p.inv_mu(dV);

    % 4-6. Update value function: 
    % Vj^(n+1) = [(rho + 1/Delta)*I - (Sj^n*Dj^n+A)]^(-1)*[u(cj^n) + 1/Delta*Vj^n]
    
    % SD 
    SD_u = spdiags(sf(:,1).*If(:,1), 0, p.I, p.I)*Df + spdiags(sb(:,1).*Ib(:,1), 0, p.I, p.I)*Db;
    SD_e = spdiags(sf(:,2).*If(:,2), 0, p.I, p.I)*Df + spdiags(sb(:,2).*Ib(:,2), 0, p.I, p.I)*Db; 

    SD = [SD_u, sparse(p.I, p.I);
          sparse(p.I, p.I), SD_e]; % 2I*2I matrix
 
    % P = SD+A
    P = SD + A; % 2I*2I matrix

    % B = [(rho + 1/dt)*I - P] 
    B = (p.rho + 1/dt)*speye(p.I*2) - P; % 2I*2I matrix

    % b = u(c) + 1/dt*V
    
    c_stacked = c(:);
    % c_stacked = reshape(c, 2*p.I, 1);
    V_stacked = V(:);

    b = p.u(c_stacked) + 1/dt*V_stacked; % 2I*1 matrix

    % V = B\b;
    V_update = B\b; % 2I*1 matrix
    V = reshape(V_update, p.I, 2);

    % SAVE OUTPUT

    if t>1
    V_t(:,:,t-1) = V;
    end

    c_t(:,:,t) = c;
    s_t(:,:,t) = zz + r_t(t) - c;
    P_t{t} = P;

end

%% 4. SOLVE KF FORWARD GIVEN A SERIES OF INTEREST RATES AND INITIAL CONDITION

for t=1:p.It

    PT= P_t{t}';

    gg_t_stacked = reshape(gg_t(:,:,t),2*p.I,1);
    gg_t_stacked_update = (speye(2*p.I) - dt*PT)\gg_t_stacked;
    
    gg_t(:,:,t+1) = reshape(gg_t_stacked_update, p.I, 2);

    %% AGGREGATE SAVING 

    S_t(t) = gg_t_stacked'*aa(:)*da;
    dS_t(t) = gg_t_stacked_update'*aa(:)*da - gg_t_stacked'*aa(:)*da;

end

%% UPDATE INTEREST RATE

    r_t_update = r_t - xi'.*dS_t;
    r_t = r_t_update;

    dist(nr) = max(abs(dS_t));

    if dist(nr)<p.tol_S
        disp('Equilibrium found. Iteration')
        disp(nr)
        break
    end
    
end

%% 6. GRAPH

%% 6-1. Initial and terminal stationary distribution

set(gca, 'FontSize', 14)
plot(a, ss_initial.gg(:,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, ss_initial.gg(:,2), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(a, ss_terminal.gg(:,1), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold on
plot(a, ss_terminal.gg(:,2), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Densities, g_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([p.amin, p.amin], yy, '--k', 'LineWidth', 2)
hold off
text(-0.1, 0.2, '$a=\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.2 1])
ylim([0 10])
legend(sprintf('Unemployed, r=%.4f', ss_initial.r), ...
       sprintf('Employed, r=%.4f', ss_initial.r), ...
       sprintf('Unemployed, r=%.4f', ss_terminal.r), ...
       sprintf('Employed, r=%.4f', ss_terminal.r), 'Location', 'best', 'FontSize', 14)
title('Initial Steady State (Solid) and Terminal Steady State (Dashed)', 'FontSize', 18)

%% 6-2. Equilibrium interest rates path

N1 = 4;
T1 = -N1*dt;
time1 = T1 + (1:N1)'*dt;
time2 = [time1;time];
r_t2 = [ss_initial.r*ones(N1,1);r_t];
set(gca,'FontSize',14)
plot(time2, r_t2, 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(time2, r_T*ones(1,N1+p.It), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold off
xlim([T1, p.tmax])
xlabel('Time')
title('Equilibrium Interest Rate, r(t)')

%% 6-3. Dynamics of wealth distribution

set(gca, 'FontSize', 14)
plot(a, ss_initial.gg(:,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, ss_initial.gg(:,2), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(a, gg_t(:, 1, 26), 'LineWidth', 2, 'LineStyle', "--", 'Color', 'r')
hold on
plot(a, gg_t(:, 2, 26), 'LineWidth', 2, 'LineStyle', "--", 'Color', 'b')
hold on
plot(a, gg_t(:, 1, 51), 'LineWidth', 2, 'LineStyle', "-.", 'Color', 'r')
hold on
plot(a, gg_t(:, 2, 51), 'LineWidth', 2, 'LineStyle', "-.", 'Color', 'b')
hold on
plot(a, ss_terminal.gg(:,1), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold on
plot(a, ss_terminal.gg(:,2), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Densities, g_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([p.amin, p.amin], yy, '--k', 'LineWidth', 2)
hold off
text(-0.11, 0.1, '$a=\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.2 0.5])
ylim([0 3])
legend({'T=0: Unemployed', 'T=0: Employed', ...
        'T=5: Unemployed', 'T=5: Employed', ...
        'T=10: Unemployed', 'T=10: Employed', ...
        'T=\infty: Unemployed', 'T=\infty: Employed'}, ...
        'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex')
title('Dynamics of Wealth Distribution', 'FontSize', 18)

%% 6-4. Dynamics of savings

set(gca, 'FontSize', 14)
plot(a, ss_initial.s(:, 1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, ss_initial.s(:, 2), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(a, s_t(:, 1, 1), 'LineWidth', 2, 'LineStyle', "--", 'Color', 'r')
hold on
plot(a, s_t(:, 2, 1), 'LineWidth', 2, 'LineStyle', "--", 'Color', 'b')
hold on
plot(a, s_t(:, 1, 6), 'LineWidth', 2, 'LineStyle', "-.", 'Color', 'r')
hold on
plot(a, s_t(:, 2, 6), 'LineWidth', 2, 'LineStyle', "-.", 'Color', 'b')
hold on
plot(a, ss_terminal.s(:, 1), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold on
plot(a, ss_terminal.s(:, 2), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Saving, s_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([p.amin, p.amin], yy, '--k', 'LineWidth', 2)
hold off
text(-0.12, -0.05, '$a=\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.2 0.5])
ylim([-0.2 0.1])
legend({'T=0: Unemployed', 'T=0: Employed', ...
        'T=0.2: Unemployed', 'T=0.2: Employed', ...
        'T=1: Unemployed', 'T=1: Employed', ...
        'T=\infty: Unemployed', 'T=\infty: Employed'}, ...
        'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex')
title('Dynamics of Saving', 'FontSize', 18)

%% 6-5. Dynamics of consumption

set(gca, 'FontSize', 14)
plot(a, ss_initial.c(:, 1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, ss_initial.c(:, 2), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(a, c_t(:, 1, 1), 'LineWidth', 2, 'LineStyle', "--", 'Color', 'r')
hold on
plot(a, c_t(:, 2, 1), 'LineWidth', 2, 'LineStyle', "--", 'Color', 'b')
hold on
plot(a, c_t(:, 1, 6), 'LineWidth', 2, 'LineStyle', "-.", 'Color', 'r')
hold on
plot(a, c_t(:, 2, 6), 'LineWidth', 2, 'LineStyle', "-.", 'Color', 'b')
hold on
plot(a, ss_terminal.c(:, 1), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold on
plot(a, ss_terminal.c(:, 2), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Consumption, c_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([p.amin, p.amin], yy, '--k', 'LineWidth', 2)
hold off
text(-0.11, 0.17, '$a=\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.2 0.5])
ylim([0.1 0.25])
legend({'T=0: Unemployed', 'T=0: Employed', ...
        'T=0.2: Unemployed', 'T=0.2: Employed', ...
        'T=1: Unemployed', 'T=1: Employed', ...
        'T=\infty: Unemployed', 'T=\infty: Employed'}, ...
        'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex')
title('Dynamics of Consumption', 'FontSize', 18)