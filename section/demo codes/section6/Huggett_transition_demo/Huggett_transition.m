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