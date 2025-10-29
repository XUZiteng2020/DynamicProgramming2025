function p = define_parameters()

% This function defines the parameters needed for the ramsey_sequential.m script

%% Economic Parameters

    % Discount rate
    p.rho = 0.03;

    % Inverse of IES 
    p.theta = 1;

    % Technology growth rate
    p.g = 0.02;

    % Population growth rate
    p.n = 0.02;

    % Capital share
    p.alpha = 1/3;

    % TFP
    p.A = 1;

%% Economic Functions
    
    % Production function
    ;

    % MPK
    ;  

%% Boundary Conditions

    % Initial capital
    p.k0 = 10;

%% Grid Paramters

    p.tmin = 0;
    p.tmax = 100;

    % The number of time steps
    p.I = 300; 

%% Newton Method Tuning Parameters
    
    % Tolerance for Newton's method
    p.tol = 1e-6;  

    % Maximum iterations for Newton's method
    % p.maxit = 1000;  

end