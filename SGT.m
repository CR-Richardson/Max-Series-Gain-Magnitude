function [alpha,data,dec]=SGT(syst)

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
% Compute the maximum series gain (alpha) when using the Small Gain Theorem (based on Lemma 6.5, Nonlinear systems, Khalil).
% Implementation assumes D = 0.
%
% Parameters:
% syst: Structure containing the system matrices of an example.
%
% Returns:
% alpha: Maximum series gain (float)
% data:  Structure containing solutions of the LMI parametrised by alpha
% dec:   # number of decision variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
A     = syst.a;
B     = syst.b;
C     = syst.c;
D     = syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

%% Initialising alpha

if m == 1
   Gm = margin(ss(A,B,-C,-D));
   if Gm > 10000
      Gm = 10000;
   end
else
    Gm = 1000;
end

% Determine initial upper/lower bound and initial test value
alpha_up  = Gm*0.999;
alpha_low = 0; % We know alpha = 0 is always feasible as system's are stable
alpha     = alpha_up;

eps       = 1e-6;
%%
% Determine alpha by repeatedly solving LMI until the largest alpha is 
% found where LMI is feasible 

while ((alpha_up - alpha_low)/alpha_up) > eps
    
    % Perform loopshift
    syst_ls = LoopShift1(syst, alpha);
    A     = syst_ls.a;
    B     = syst_ls.b;
    C     = syst_ls.c;
    D     = syst_ls.d;

    setlmis([]);
    
    P   = lmivar(1,[n,1]);
    k = 1.0; % sector bound
    
    % LMI
    lmiterm([1,1,1,P],A',1,'s');
    lmiterm([1,1,1,0],k*(C'*C) ); 
    
    lmiterm([1,1,2,P],1,B);
    
    lmiterm([1,2,2,0],-k);
    
    % P > 0
    lmiterm([2,1,1,P],-1,1);
    
    LMISYS=getlmis;
    [tmin,xfeas] = feasp(LMISYS,[1e-20 5000 -0.1 1000 1]);
    
    % Update alpha upper/lower bound plus new test value 
     if tmin < 0 % if LMIs are feasible
        alpha_low = alpha;
     else 
        alpha_up = alpha; % if LMIs are infeasible
     end
      
    alpha = (alpha_up + alpha_low)/2;
    
    if alpha<1e-9
        break;
    end
end

%% Return solutions

if alpha<1e-9
    disp('Solution cannot be found!');
    dec       =  decnbr(LMISYS); % returns number of decision varibles
    data.P    =  nan;
    alpha     =  nan;
else
    dec       =  decnbr(LMISYS); % returns number of decision varibles
    data.P    =  dec2mat(LMISYS,xfeas,P);
end

end


