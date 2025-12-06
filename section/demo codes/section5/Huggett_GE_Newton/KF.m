function [gg] = KF(r, p, G, P);

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
    g_sum = g_stacked'*ones(2*p.I,1)*G.da;
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


end