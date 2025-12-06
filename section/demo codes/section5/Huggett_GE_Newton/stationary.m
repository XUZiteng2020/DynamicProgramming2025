function [S, ss] = stationary(r, p, G)

%% HJB

[c, s, V, P] = HJB(r, p, G);

%% KF

[gg] = KF(r, p, G, P);

%% aggregate saving S(r)
    
S = gg(:,1)'*G.a*G.da + gg(:,2)'*G.a*G.da;

%% OUTPUT

ss.c = c;
ss.s = s;
ss.V = V;
ss.gg = gg;

end