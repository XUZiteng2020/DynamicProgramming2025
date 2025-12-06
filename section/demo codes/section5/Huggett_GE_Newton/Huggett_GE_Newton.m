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

G.a = linspace(p.amin, p.amax, p.I)';
G.da = (p.amax-p.amin)/(p.I-1);

G.aa = [G.a, G.a]; % I*2 matrix

%% NEWTON'S METHOD

% Guess an initial value of interest rate
r0 = 0.03;

% Define aggregate saving S(r) as a function of r
S = @(r) stationary(r, p, G);

% Solve for equilibrium interest rate r such that S(r)=0 using fsolve
options = optimoptions("fsolve", "FunctionTolerance", p.tol_S);
r_equilibrium = fsolve(S, r0, options);

% Solve for stationary equilibrium given the equilibrium interest rate
[~, ss] = stationary(r_equilibrium, p, G);

%% 6. GRAPHS 

%% 6-1. Optimal consumption 

set(gca, 'FontSize', 18)
plot(G.a, ss.c, 'LineWidth', 2)
grid
xlabel('Wealth, a','FontSize', 14)
ylabel('Consumption, c_j(a)','FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

%% 6-2. Optimal savings 

set(gca, 'FontSize', 18)
plot(G.a, ss.s, G.a, zeros(1,p.I), '--k', 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Saving, s_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

%% 6-3. Value function

set(gca, 'FontSize', 18)
plot(G.a, ss.V, 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Value function, V_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

%% 6-4. Wealth distribution

set(gca, 'FontSize', 14)
plot(G.a, ss.gg, 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Densities, g_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([-0.15,-0.15], yy, '--k', 'LineWidth', 2)
text(-0.15, yy(1)-0.02*(yy(2) - yy(1)), '$\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.2 1])
legend('Unemployed', 'Employed', 'Borrowing Constraint', 'Location', 'best', 'FontSize', 14)