function  [c, s, V, P] = HJB(r, p, G);

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the forward and backward finite-difference operators  
% Df such that Df*V=dVf and Db such that Db*V=dVb

    Df = zeros(p.I, p.I);
    Db = zeros(p.I, p.I);

    for i=1:p.I-1
        Df(i,i) = -1/G.da; Df(i,i+1) = 1/G.da;
    end
    Df = sparse(Df);

    for i=2:p.I
        Db(i,i-1) = -1/G.da; Db(i,i) = 1/G.da;
    end
    Db = sparse(Db);

% 3-2. Construct A matrix

    A = [-p.lambda(1).*speye(p.I), p.lambda(1).*speye(p.I);
         p.lambda(2)*speye(p.I)  , -p.lambda(2)*speye(p.I)]; % 2I*2I matrix

% 3-3. Guess an initial value of the value function
 
    % The value function of "staying put" 
    
    zz = ones(p.I,1).*p.zz; % I*2 matrix
    V0 = p.u(zz + r.*G.aa)/p.rho; % I*2 matrix; first column is for unemployed; second clumn is for employed;
    V = V0;

%% 4. VALUE FUNCTION ITERATION


for n=1:p.maxit

    % 4-1. Compute the derivative of the value function 
    dVf = Df*V; % I*2 matrix
    dVb = Db*V;

    % 4-2. Boundary conditions
    % a>=a_min is enforced (borrowing constraint)
    dVb(1, :) = p.mu(p.zz + r*G.aa(1,:));
   
    % a<=a_max is enforced which helps stability of the algorithm
    dVf(end, :) = p.mu(p.zz + r*G.aa(end,:));

    I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)

    % 4-3. Compute the optimal consumption
    cf = p.inv_mu(dVf); % I*2 matrix
    cb = p.inv_mu(dVb);
    
    % 4-4. Compute the optimal savings
    sf = zz + r.*G.aa - cf; % I*2 matrix
    sb = zz + r.*G.aa - cb;

    % 4-5. Upwind scheme
    If = sf>0;  % I*2 matrix
    Ib = sb<0;
    dV0 = p.mu(zz + r.*G.aa);

    dV = dVf.*If + dVb.*Ib + dV0.*(1-If-Ib);  % I*2 matrix

    c = p.inv_mu(dV);
    s = zz + r.*G.aa - c;

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

end
