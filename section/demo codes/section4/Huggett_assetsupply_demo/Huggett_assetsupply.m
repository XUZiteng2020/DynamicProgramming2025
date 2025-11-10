%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_Huggett_assetsupply
% 
% Author: Kiyea Jin
% Date: Nov 8, 2025
%
% Description:
% This MATLAB script solves the HJB equation and the KF equation 
% of the Huggett model. Also, it computes the asset supply S(r) for a
% series of r.
%
% Reference: Huggett_partialeq.m by Benjamin Moll
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-sigma))/(1-sigma)
% - Relative risk aversion coefficient (sigma): 2
% - Discrete grid of interest rates levels (rgrid): -0.05 to 0.04
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

p = define_parameters_assetsupply();






%% 6. GRAPHS 

set(gca,'FontSize',14)
plot(dist,'LineWidth',2)
grid
xlabel('Iteration')
ylabel('||V^{n+1} - V^n||')

% Verr = c.^(1-s)/(1-s) + dV_Upwind.*(zz + r.*aa - c) + ones(I,1)*la.*(V_switch - V) - rho.*V;
% 
% set(gca,'FontSize',14)
% plot(a,Verr,'LineWidth',2)
% grid
% xlabel('k')
% ylabel('Error in HJB Equation')
% xlim([amin amax])

% 6-1. Optimal consumption 

set(gca, 'FontSize', 18)
plot(a, c_r(:,1,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, c_r(:,2,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(a, c_r(:,1,19), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold on
plot(a, c_r(:,2,19), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b')
hold off
grid
xlabel('Wealth, a','FontSize', 14)
ylabel('Consumption, c_j(a)','FontSize', 14)
xlim([p.amin p.amax])
legend(sprintf('Unemployed, r=%.4f', rgrid(1)), ...
       sprintf('Employed, r=%.4f', rgrid(1)), ...
       sprintf('Unemployed, r=%.4f', rgrid(19)), ...
       sprintf('Employed, r=%.4f', rgrid(19)), 'Location', 'best', 'FontSize', 14)

% 6-2. Optimal savings 

set(gca, 'FontSize', 18)
plot(a, adot_r(:,1,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, adot_r(:,2,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(a, adot_r(:,1,19), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold on
plot(a, adot_r(:,2,19), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Saving, s_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend(sprintf('Unemployed, r=%.4f', rgrid(1)), ...
       sprintf('Employed, r=%.4f', rgrid(1)), ...
       sprintf('Unemployed, r=%.4f', rgrid(19)), ...
       sprintf('Employed, r=%.4f', rgrid(19)), 'Location', 'best', 'FontSize', 14)

% 6-3. Value function

set(gca, 'FontSize', 18)
plot(a, V_r(:,1,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, V_r(:,2,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(a, V_r(:,1,19), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold on
plot(a, V_r(:,2,19), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b')
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Value function, V_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend(sprintf('Unemployed, r=%.4f', rgrid(1)), ...
       sprintf('Employed, r=%.4f', rgrid(1)), ...
       sprintf('Unemployed, r=%.4f', rgrid(19)), ...
       sprintf('Employed, r=%.4f', rgrid(19)), 'Location', 'best', 'FontSize', 14)

% 6-4. Wealth distribution

set(gca, 'FontSize', 14)
plot(a, g_r(:,1,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(a, g_r(:,2,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
plot(a, g_r(:,1,19), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r')
hold on
plot(a, g_r(:,2,19), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b')
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
legend(sprintf('Unemployed, r=%.4f', rgrid(1)), ...
       sprintf('Employed, r=%.4f', rgrid(1)), ...
       sprintf('Unemployed, r=%.4f', rgrid(19)), ...
       sprintf('Employed, r=%.4f', rgrid(19)), 'Location', 'best', 'FontSize', 14)

%% 6-5. ASSET SUPPLY FUNCTION S(R)

set(gca,'FontSize',14)
plot(S, rgrid, 'LineWidth', 2, 'Color', 'r')
hold on;
aaa = linspace(p.amin-0.01, p.amax, p.Ir)';
plot(aaa, p.rho*ones(p.Ir, 1), 'LineWidth', 2, 'Linestyle', '--', 'Color', 'k')
hold on;
rrr = linspace(p.rmin, p.rho+0.01, p.Ir)';
plot(p.amin*ones(p.Ir, 1), rrr, 'LineWidth', 2, 'Linestyle', '--', 'Color', 'k')
hold on;
plot(zeros(p.Ir, 1), rrr, 'LineWidth', 2, 'Color', 'b');
ylabel('Interest rate, r','FontSize',16)
xlabel('Asset supply, S(r)','FontSize',16)
ylim([p.rmin p.rho+0.01])
maxS = max(S);
xlim([p.amin-0.01 maxS])
text(-0.1,0.045,'$r = \rho$','FontSize',16,'interpreter','latex')
text(-0.07,0,'$S(r)$','FontSize',16,'interpreter','latex')
text(-0.145,0.01,'$a=\underline{a}$','FontSize',16,'interpreter','latex')
text(0.005,0,'$B$','FontSize',16,'interpreter','latex')
