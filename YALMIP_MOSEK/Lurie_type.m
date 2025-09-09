%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% CR Richardson and MC Turner
% ECS
% University of Southampton
% UK
%
% Date: 5/09/25
%
% Purpose: 
% Compute the maximum series gain (alpha) when using the Lurie type Lyapunov
% function as defined by Corollary 2.
%
% Parameters:
% syst:      Structure containing the system matrices of an example.
% eps:       Loop termination accuracy.
% alpha_up:  Largest series gain to be tested.
% alpha_low: Lowest series gain to be tested.
%
% Returns:
% alpha_low: Maximum series gain.
% data:      Structure containing solutions of the LMI parametrised by alpha.
% dec:       # of decision variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha_low,data,dec]=Lurie_type(syst,eps,alpha_up,alpha_low)

%% Parameters
B     = syst.b;
D     = syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

% Raise error if D != 0.
if ~isequal(D, zeros(m))
    error('D terms must be added to LMIs.');
end

%% Initialising alpha
alpha = alpha_up;

%% Determine max alpha by repeatedly solving LMI for largest feasible alpha 

while ((alpha_up-alpha_low)/alpha_up) > eps

% Perform loopshift
syst_ls = LoopShift1(syst, alpha);
A     = syst_ls.a;
B     = syst_ls.b;
C     = syst_ls.c;
D     = syst_ls.d;

% define variables
P = sdpvar(n,n,'symmetric');
L1 = sdpvar(m,m,'diagonal');
L2 = sdpvar(m,m,'diagonal');
V = sdpvar(m,m,'full');
Q13 = sdpvar(m,m,'full');
Q22 = sdpvar(m,m,'symmetric');

% define constraints
epsilon = 1e-2;
In = eye(n);
Im = eye(m);
c1 = P >= epsilon*In; % PD
c2 = L1 >= epsilon*Im; % PD
c3 = L2 >= epsilon*Im; % PD
c4 = V(~Im) >= 0; % Metzler
c5 = Q13(:) >= 0; % Non-negative matrix
c6 = Q22(:) >= 0; % Non-negative matrix

X11 = ( P*A + C'*(L1 + L2)*C*A - C'*(V - Q22)*C ) + ( P*A + C'*(L1 + L2)*C*A - C'*(V - Q22)*C )';
X21 = ( B'*P + B'*C'*(L1 + L2)*C - V*C + V'*C + (2*Q22 - Q13)'*C + (L1 - L2)*C*A );
X22 = ( V + Q22 + Q13 + (L1 - L2)*C*B ) + ( V + Q22 + Q13 + (L1 - L2)*C*B )';

c7 = [X11, X21';
      X21, X22] <= -epsilon*eye(n+m); % ND

Constraints = [c1, c2, c3, c4, c5, c6, c7];

% solve LMI
Objective = [];
options = sdpsettings('solver', 'mosek', 'verbose', 0);
result = optimize(Constraints, Objective, options);

% Update alpha upper/lower bound plus new test value
if result.problem == 0 % if LMIs are feasible
    alpha_low = alpha;
else 
    alpha_up = alpha; % if LMIs are infeasible
end
  
alpha = (alpha_up + alpha_low)/2;
end

%% Output

% Get variable indices
vars = [P(:); L1(:); L2(:); V(:); Q13(:); Q22(:)];
dec = length(unique(getvariables(vars)));

% solutions
data.P = value(P);
data.L1 = value(L1);
data.L2 = value(L2);
data.V = value(V);
data.Q13 = value(Q13);
data.Q22 = value(Q22);

end