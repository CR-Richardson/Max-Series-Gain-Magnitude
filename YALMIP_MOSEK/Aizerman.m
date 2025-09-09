
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% MC Turner and CR Richardson
% ECS
% University of Southampton
% UK
%
% Date: 20/05/25
%
% Purpose: 
% Computes the Nyquist gain (alpha) for the example systems
% based on the Aizerman conjecture. The algorithm assumes the set of
% stabilising gains is connected, which is often, but not always the case.
%
% Scripts
% Examples: Contains example linear systems.
%
% Variables:
% Total_Ex:    Total number of examples.
% Ex_array:    The set of example systems to compute alpha for.
% eps:         Loop termination accuracy.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script variables
clear all; close all;
Total_Ex   = 12;
Ex_array   = 9:12;
eps        = 1e-6;

%% Makes example systems accessible to script
Examples;

%% Find largest slope size for each example

alpha = zeros(1,Total_Ex);

for i=Ex_array

    % calculate initial upper and lower slope size
    curlyA_u=100000;         % Upper value of slope size
    curlyA_l=0;             % Lower value of slope size (we know 0 is always feasible as system is stable)
    curlyA=curlyA_u;

    % Load system
    A     = Syst{i}.a;
    B     = Syst{i}.b;
    C     = Syst{i}.c;
    D     = Syst{i}.d;
    [n,m] = size(B);

    % Raise error if D != 0. 
    if ~isequal(D, zeros(m))
        error('D is not a zero matrix of size %d x %d.', m, m);
    end
    
    % Determine largest slope size (curlyA)
    while ((curlyA_u-curlyA_l)/curlyA_u) > eps
        
        Acl = A + curlyA*B*C;
        E   = eig(Acl);
        rho = max(real(E));

        if rho < 0
            curlyA_l = curlyA;        % System is stable
        else 
            curlyA_u = curlyA;        % System is unstable
        end
  
        curlyA = (curlyA_u + curlyA_l) / 2;

    end

    alpha(i) = curlyA_l;

end

for i=Ex_array
    disp(['Example ', num2str(i), ' ', num2str(alpha(i))]);
end

