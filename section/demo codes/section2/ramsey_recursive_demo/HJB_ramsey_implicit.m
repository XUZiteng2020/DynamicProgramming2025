%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_ramsey_implicit
% 
% Author: Kiyea Jin
% Date: Oct 28, 2025
%
% Description:
% This MATLAB script implements implicit method to solve the HJB equation
% of the deterministic Neoclassical Growth Model using mixed method. 
% Boundary conditions are not yet considered in this script.
%
% Reference:
% HJB_NGM_implicit.m by Benjamin Moll
% ramsey_implicit.m by Pontus Rendahl
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-gamma))/(1-gamma)
% - Production function: f(k) = A*k^alpha
% - Relative risk aversion coefficient (gamma): 2
% - Discount rate (rho): 0.03
% - Depreciation rate (delta): 0.025
% - Elasticity of output with respect to capital (alpha): 1/3
% - Total fator productivity (A): 1
% - Delta = 1000 (Can be arbitrarily large in implicit method)
%
% Code Structure:
% 1. DEFINE PARAMETERS
% 2. INITIALIZE GRID POINTS
% 3. PRE-ITERATION INITIALIZATION
% 4. VALUE FUNCTION ITERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters();

%% 2. INITIALIZE GRID POINTS

% Steady-state level of capital: f'(kss)=rho+delta
;

% The parameter klim controls the log-distance from kss
;
;

% Construct grid points
;
;

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the differential operator D such that DV=dV

    D = zeros(p.I,p.I);
    
    % Forward differencing for i=1
    ;
    ;
    
    % Backward differencing for i=I
    ;
    ;
    
    % Central differencing for i=2,...,I-1
    for 
        ; 
        ;
    end

% 3-2. Guess an initial value of the value function

    V0 = p.u(p.f(k)) / p.rho;
    V = V0;    

%% 4. VALUE FUNCTION ITERATION

tic;

for 

    % 4-1. Compute the derivative of the value function
    ;
        
    % 4-2. Compute the optimal consumption
    ;

    % 4-3. Compute the optimal savings
    ;

    % 4-4. Update the value function: V^(n+1) = [(rho+1/Delta)*I - SD]^(-1)[u(c) + 1/Delta*V^n]
    
        % B = [(rho+1/Delta)*I - SD]
        ;
        ;

        % b = [u(c) + 1/Delta*V^n]
        ;

        % V^(n+1) = B\b
        ;
        
        % Update the value function
        ;
        ;


    % 4-5. Check convergence
          
        ;

        if 
        disp('Value function converged. Iteration = ')
        disp(n)
        break
        end
end

toc;

%% 5. Saddle path

% Saddle path
figure;
p1 = plot(k, c, 'linewidth', 2);
set(gca, 'FontSize', 18)
xlabel('Capital, k','FontSize', 18)
ylabel('Consumption, c','FontSize',18)

% dk/dt = 0
k_null_c = p.f(k) - p.delta*k;
hold on
p2 = plot(k, k_null_c, 'linewidth', 2);
yy = get(gca, 'yLim');

% dc/dt = 0
p3 = plot([kss kss], yy, 'linewidth', 2);

legend1 = legend([p1,p2,p3], 'Saddle path', '\Delta k=0', '\Delta c=0');
set(legend1, 'Location', 'best', 'FontSize', 18)
hold off
