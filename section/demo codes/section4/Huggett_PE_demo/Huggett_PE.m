%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_Huggett_PE
% 
% Author: Kiyea Jin
% Date: Nov 8, 2025
%
% Description:
% This MATLAB script solves the HJB equation and the KF equation 
% of the Huggett model.
%
% Reference: Huggett_partialeq.m by Benjamin Moll
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-sigma))/(1-sigma)
% - Relative risk aversion coefficient (sigma): 1.2
% - Interest rate (r) : 0.035
% - Discount rate (rho): 0.05
% - Income: z = [z_u, z_e] = [0.1, 0.2];
% - Lambda: la = [la_u, la_e] = [1.5, 1];
% - Discrete grid of asset levels (a): -0.02 to 2
% - Borrowing constraint: a>=-0.02
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

p = define_parameters();

%% 2. INITIALIZE GRID POINTS


%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the forward and backward finite-difference operators  
% Df such that Df*V=dVf and Db such that Db*V=dVb


% 3-2. Construct A matrix


% 3-3. Guess an initial value of the value function
 
    % The value function of "staying put" 
    

%% 4. VALUE FUNCTION ITERATION

for n=1:p.maxit

    % 4-1. Compute the derivative of the value function 


    % 4-2. Boundary conditions
    % a>=a_min is enforced (borrowing constraint)
   
    % a<=a_max is enforced which helps stability of the algorithm

    I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)

    % 4-3. Compute the optimal consumption

    
    % 4-4. Compute the optimal savings


    % 4-5. Upwind scheme


    % 4-6. Update value function: 
    % Vj^(n+1) = [(rho + 1/Delta)*I - (Sj^n*Dj^n+A)]^(-1)*[u(cj^n) + 1/Delta*Vj^n]
    
    % SD 

 
    % P = SD+A


    % B = [(rho + 1/Delta)*I - P] 


    % b = u(c) + 1/Delta*V
    

    % V = B\b;


    % 3-6. Convergence criterion


    if 
       disp('Value function converged. Iteration = ')
       disp(n)
       break
    end
end

toc;

%% 5. KF EQUATION

%% 5-1. SOLVE FOR 0=gdot=P'*g

    % Define P'

    % Define gdot

    % Since the eigenvector is only defined up to a scalar, we need to fix one
    % value; otherwise matrix is singular.

    
    % g=P'\gdot

    % Normalization


    % Reshape

    
%% 5-2. KF EQUATION SOLVED WITH EIGS
% Notes: [V, D] = eigs(A, k, sigma) returns the k largest (in magnitude) eigenvalues D and 
% corresponding eigenvectors V of matrix A closest to the target sigma. 
% - A: The input matrix.
% - k: The number of eigenvalues (and corresponding eigenvectors) to compute.
% - If SIGMA is a real or complex scalar including 0, eigs finds the eigenvalues closest to SIGMA.
    
    PT_eigs = P';
    [g_stacked_eigs, eigval] = eigs(PT_eigs, 1, 0);
    g_sum_eigs = g_stacked_eigs'*ones(2*p.I,1)*da;
    g_stacked_eigs = g_stacked_eigs./g_sum_eigs;

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

%% 6-1. Optimal consumption 

set(gca, 'FontSize', 18)
plot(a, c, 'LineWidth', 2)
grid
xlabel('Wealth, a','FontSize', 14)
ylabel('Consumption, c_j(a)','FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

%% 6-2. Optimal savings 

adot = zz + p.r.*aa - c;

set(gca, 'FontSize', 18)
plot(a, adot, a, zeros(1,p.I), '--k', 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Saving, s_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

%% 6-3. Value function

set(gca, 'FontSize', 18)
plot(a, V, 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Value function, V_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

%% 6-4. Wealth distribution

set(gca, 'FontSize', 14)
plot(a, gg, 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Densities, g_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([-0.02,-0.02], yy, '--k', 'LineWidth', 2)
text(-0.02, yy(1)-0.02*(yy(2) - yy(1)), '$\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.03 1])
legend('Unemployed', 'Employed', 'Borrowing Constraint', 'Location', 'best', 'FontSize', 14)