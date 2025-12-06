%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_Huggett_GE
% 
% Author: Kiyea Jin
% Date: Nov 9, 2025
%
% Description:
% This MATLAB script solves the general equilibrium of the Huggett model,
% finding the equilibrium interest rate that clears the bond market.
%
% Reference: Huggett_equilibrium_iterate.m by Benjamin Moll
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-sigma))/(1-sigma)
% - Elasticity of intertemporal substitution (sigma): 2
% - Discount rate (rho): 0.05
% - Income: z = [z_u, z_e] = [0.1, 0.2];
% - Lambda: la = [la_u, la_e] = [1.2, 1.2];
% - Discrete grid of asset levels (a): -0.15 to 5
% - Borrowing constraint: a>=-0.15
% - Delta = 1000; (Can be arbitrarily large in implicit method)
%
% Code Structure:
% 1. DEFINE PARAMETERS
% 2. INITIALIZE GRID POINTS
% 3. PRE-ITERATION INITIALIZATION
% 4. VALUE FUNCTION ITERATION
% 5. KF EQUATION
% 6. GRAPHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters_GE();

%% 2. INITIALIZE GRID POINTS

a = linspace(p.amin, p.amax, p.I)';
da = (p.amax-p.amin)/(p.I-1);

aa = [a, a]; % I*2 matrix

% %% INITIALIZE GRID POINTS FOR INTEREST RATES
% 
% rgrid = linspace(p.rmin, p.rmax, p.Ir)'; % Ir*1 matrix

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

% GUESS AN INITIAL VALUE OF THE INTEREST RATE

    r0 = 0.03;
    rmin0 = 0.01;
    rmax0 = 0.04;

    r = r0;
    r_min = rmin0;
    r_max = rmax0;

% 3-3. Guess an initial value of the value function
 
    % The value function of "staying put" 
    
    zz = ones(p.I,1).*p.zz; % I*2 matrix
    V0 = p.u(zz + r0.*aa)/p.rho; % I*2 matrix; first column is for unemployed; second clumn is for employed;
    V = V0;

%% 4. VALUE FUNCTION ITERATION

for nr=1:p.Nr

    r_r(nr) = r;
    rmin_r(nr) = r_min;
    rmax_r(nr) = r_max;

    % Use the value function solution from the previous interest rate iteration 
    % as the initial guess for the next iteration

    if nr>1
       V0 = V_r(:,:,nr-1);
       V = V0;
    end

for n=1:p.maxit

    % 4-1. Compute the derivative of the value function 
    dVf = Df*V; % I*2 matrix
    dVb = Db*V;

    % 4-2. Boundary conditions
    % a>=a_min is enforced (borrowing constraint)
    dVb(1, :) = p.mu(p.zz + r*aa(1,:));
   
    % a<=a_max is enforced which helps stability of the algorithm
    dVf(end, :) = p.mu(p.zz + r*aa(end,:));

    I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)

    % 4-3. Compute the optimal consumption
    cf = p.inv_mu(dVf); % I*2 matrix
    cb = p.inv_mu(dVb);
    
    % 4-4. Compute the optimal savings
    sf = zz + r.*aa - cf; % I*2 matrix
    sb = zz + r.*aa - cb;

    % 4-5. Upwind scheme
    If = sf>0;  % I*2 matrix
    Ib = sb<0;
    dV0 = p.mu(zz + r.*aa);

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

    % B = [(rho + 1/Delta)*I - P] 
    B = (p.rho + 1/p.Delta)*speye(p.I*2) - P; % 2I*2I matrix

    % b = u(c) + 1/Delta*V
    
    c_stacked = c(:);
    % c_stacked = reshape(c, 2*p.I, 1);
    V_stacked = V(:);

    b = p.u(c_stacked) + 1/p.Delta*V_stacked; % 2I*1 matrix

    % V = B\b;
    V_update = B\b; % 2I*1 matrix

    V_change = V_update - V_stacked;
    V = reshape(V_update, p.I, 2);

    % 3-6. Convergence criterion
    dist(n) = max(abs(V_change));

    if dist(n)<p.tol
       disp('Value function converged. Iteration = ')
       disp(n)
       break
    end
end

toc;

%% 5. KF EQUATION

%% 5-1. SOLVE FOR 0=gdot=P'*g

    % Define P'
    PT = P';

    % Define gdot
    gdot_stacked = zeros(2*p.I, 1);

    % Since the eigenvector is only defined up to a scalar, we need to fix one
    % value; otherwise matrix is singular.
    ifix = 1;
    gdot_stacked(ifix) = 0.1;

    rowfix = [zeros(1,ifix-1), 1, zeros(1, 2*p.I-ifix)];
    PT(ifix, :) = rowfix;
    
    % g=P'\gdot
    g_stacked = PT\gdot_stacked; % 2I*1 matrix

    % Normalization
    g_sum = g_stacked'*ones(2*p.I,1)*da;
    g_stacked = g_stacked./g_sum;

    % Reshape
    gg = reshape(g_stacked, p.I, 2);

% %% 5-2. KF EQUATION SOLVED WITH EIGS
% % Notes: [V, D] = eigs(A, k, sigma) returns the k largest (in magnitude) eigenvalues D and 
% % corresponding eigenvectors V of matrix A closest to the target sigma. 
% % - A: The input matrix.
% % - k: The number of eigenvalues (and corresponding eigenvectors) to compute.
% % - If SIGMA is a real or complex scalar including 0, eigs finds the eigenvalues closest to SIGMA.
% 
%     PT_eigs = P';
%     [g_stacked_eigs, eigval] = eigs(PT_eigs, 1, 0);
%     g_sum_eigs = g_stacked_eigs'*ones(2*p.I,1)*da;
%     g_stacked_eigs = g_stacked_eigs./g_sum_eigs;

%% COMPUTE VARIABLES FOR A GIVEN r(ir) for all ir=1:p.Ir
% Notes: Each matrix has dimensions p.I*2(u,e)*p.Ir
    
    % c, adot, V, g
    c_r(:,:,nr) = c;
    adot_r(:,:,nr) = zz + r.*aa - c;
    V_r(:,:,nr) = V;
    g_r(:,:,nr) = gg;
    
    % aggregate saving S(r)
    S(nr) = gg(:,1)'*a*da + gg(:,2)'*a*da;

%% UPDATE INTEREST RATE

    if S(nr)>p.tol_S
       disp('Excess supply')
       % Decrease r if S(r)>0
       r_max = r;
       r = (r_min+r_max)/2;           

    elseif S(nr)<-p.tol_S
       disp('Excess demand')
       % Increase r if S(r)<0
       r_min = r;
       r = (r_min+r_max)/2;
        
    elseif abs(S(nr))<p.tol_S
       disp('Equilibrium found. Iteration = ')
       disp(n)
       break
    end
   
end

%% 6. GRAPHS 

% 6-1. Optimal consumption 

set(gca, 'FontSize', 18)
plot(a, c_r(:,1,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, c_r(:,2,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold off
grid
xlabel('Wealth, a','FontSize', 14)
ylabel('Consumption, c_j(a)','FontSize', 14)
xlim([p.amin p.amax])
legend(sprintf('Unemployed, r=%.4f', r), ...
       sprintf('Employed, r=%.4f', r), 'Location', 'best', 'FontSize', 14)

% 6-2. Optimal savings 

set(gca, 'FontSize', 18)
plot(a, adot_r(:,1,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, adot_r(:,2,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Saving, s_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend(sprintf('Unemployed, r=%.4f', r), ...
       sprintf('Employed, r=%.4f', r), 'Location', 'best', 'FontSize', 14)

% 6-3. Value function

set(gca, 'FontSize', 18)
plot(a, V_r(:,1,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, V_r(:,2,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Value function, V_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend(sprintf('Unemployed, r=%.4f', r), ...
       sprintf('Employed, r=%.4f', r), 'Location', 'best', 'FontSize', 14)

% 6-4. Wealth distribution

set(gca, 'FontSize', 14)
plot(a, g_r(:,1,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, g_r(:,2,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Densities, g_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([p.amin, p.amin], yy, '--k', 'LineWidth', 2)
hold off
text(-0.15, yy(1)-0.02*(yy(2) - yy(1)), '$\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.2 1])
legend(sprintf('Unemployed, r=%.4f', r), ...
       sprintf('Employed, r=%.4f', r), 'Location', 'best', 'FontSize', 14)