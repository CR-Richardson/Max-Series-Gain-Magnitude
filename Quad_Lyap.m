function [alpha,data,dec]=Quad_Lyap(syst)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% MC Turner and CR Richardson 
% ECS
% University of Southampton
% UK
%
% Date: 25/06/24
%
% Purpose: 
% Compute the maximum series gain (alpha) when using the Quadratic Criterion (Theorem 1).
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

%%
% Determine alpha by repeatedly solving LMI until the largest alpha is 
% found where LMI is feasible 

while ((alpha_up - alpha_low)/alpha_up) > 0.0001
    
    % Perform loopshift
    syst_ls = LoopShift1(syst, alpha);
    A     = syst_ls.a;
    B     = syst_ls.b;
    C     = syst_ls.c;
    D     = syst_ls.d;

    setlmis([]);
    
    P   = lmivar(1,[n,1]);
    V   = lmivar(2,[m,m]);
    Q13 = lmivar(2,[m,m]);
    Q22 = lmivar(1,[m,1]);
    
    % LMI
    lmiterm([1,1,1,P],A',1,'s');
    lmiterm([1,1,1,V],-C',C,'s');
    lmiterm([1,1,1,Q22],C',C,'s');
    
    lmiterm([1,1,2,P],1,B);
    lmiterm([1,1,2,-V],-C',1);
    lmiterm([1,1,2,V],C',1);
    lmiterm([1,1,2,Q22],2*C',1);
    lmiterm([1,1,2,Q13],-C',1);
    
    lmiterm([1,2,2,V],1,1,'s');
    lmiterm([1,2,2,Q22],1,1,'s');
    lmiterm([1,2,2,Q13],1,1,'s');
    
    % P > 0
    lmiterm([2,1,1,P],-1,1);
    
    % V: Metzler matrix conditions
    JJ = eye(m);
    count = 3;     
    for i = 1:m
        for j = 1:m
            e1 = JJ(i,:); e2 = JJ(:,j);
            if i ~= j
               lmiterm([count,1,1,V],-0.5*e1,e2,'s');
               count = count+1;
            end
        end
    end
    
    % Q13: Positivity matrix conditions
    for i = 1:m
        for j = 1:m
            e1 = JJ(i,:); e2 = JJ(:,j);
            lmiterm([count,1,1,Q13],-0.5*e1,e2,'s');
            count = count+1;
         end
    end 
    
    % Q22: Positivity matrix conditions
    for i = 1:m
        for j = 1:m
            e1 = JJ(i,:); e2 = JJ(:,j);
            lmiterm([count,1,1,Q22],-0.5*e1,e2,'s');
            count = count+1;
         end
    end 
    
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
    data.V    =  nan;
    data.Q13  =  nan;
    data.Q22  =  nan;
    alpha     =  nan;
else
    dec       =  decnbr(LMISYS); % returns number of decision varibles
    data.P    =  dec2mat(LMISYS,xfeas,P);
    data.V    =  dec2mat(LMISYS,xfeas,V);
    data.Q13  =  dec2mat(LMISYS,xfeas,Q13);
    data.Q22  =  dec2mat(LMISYS,xfeas,Q22);
end

end


