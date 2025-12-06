function [ss] = Huggett_GE(p)

%% 2. INITIALIZE GRID POINTS

G.a = linspace(p.amin, p.amax, p.I)';
G.da = (p.amax-p.amin)/(p.I-1);

G.aa = [G.a, G.a]; % I*2 matrix

%% NEWTON'S METHOD

% 1) Guess an initial value for interest rate
r0_guess = 0.03;

% 2) Define a function S(r) that computes aggregate saving
S = @(r) stationary(r, p, G);

% 3) Solve for equilibrium interest rate r such that S(r)=0 using fsolve
options = optimoptions("fsolve", "FunctionTolerance", p.tol_S);
r_equilibrium = fsolve(S, r0_guess, options);

% 4) Solve for stationary equilbrium (c, s, V, g) using r_equilibrium
[~, ss] = stationary(r_equilibrium, p, G);
ss.r = r_equilibrium;

end