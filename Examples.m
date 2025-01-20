%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% MC Turner and CR Richardson 
% ECS
% University of Southampton
% UK
%
% Date: 15/05/23
%
% Purpose: 
% Script containing the example linear systems for the max. series gain
% experiments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Park - IEEETAC 2002 - Example 3
A        = -diag([1 4 6 2 9 8 3 10 12]);
B        = -[1 0 0 0 1 0 0 1 0;
             0 1 0 1 0 0 0 0 1;
             0 0 1 0 0 1 1 0 0]';
C        =  [1 1 1 0 0 0 0 0 0;
             0 0 0 1 1 1 0 0 0;
             0 0 0 0 0 0 1 1 1];
D        =  zeros(3,3);
Syst{1}  = ss(A,B,C,D);
clear A B C D;

%% Drummond, Guiver, and Turner - IEEETAC 2022 - Example 3
A       = diag([-1 -10 -100])+0.01*ones(3,3);
B       = eye(3);
C       = [0 1 1
           0 0 1
           0 0 1];
D       = zeros(3,3);
Syst{2} = ss(A,B,C,D);
clear A B C D;

%% Fetzer and Scherer - IJRNC 2017 - Example 4.9
A       = [-4 -3  0
            2  0  0
           -1 -1 -2];
B       = [0 4 1 3
           2 0 3 1
           1 0 3 1];
C       = [-0.1 -0.2 1.0
           -1.0 -0.3 0.1
           -0.2  0.1 1.0
            0.1 -0.2 0.2];
D       = zeros(4,4);
Syst{3} = ss(A,B,C,D);
clear A B C D;

%% Turner, Kerr, and Sofrony - IJRNC 2015 - Example 22
A        = [-5.8097   -1.7258   -1.3496   -0.8882   -0.3991   -0.2843   -0.3009   -0.3524
             8.0000         0         0         0         0         0         0         0
             0    2.0000         0         0         0         0         0         0
             0         0    2.0000         0         0         0         0         0
             0         0         0    2.0000         0         0         0         0
             0         0         0         0    1.0000         0         0         0
             0         0         0         0         0    0.5000         0         0
             0         0         0         0         0         0    0.2500         0];
B        = [-0.0697    0.1982   -0.4746   -0.3289
            -0.4567    0.3788   -0.2270    0.5888
             0.3737   -0.6601    0.6961   -0.3923
            -0.7929   -0.6439   -0.0461   -0.2640
             0.0148    0.4926    0.5402    0.7878
            -0.5594   -0.0041    0.2616    0.4186
             0.0385   -0.5487   -0.5681    0.5668
            -0.1095    0.0133    0.6951    0.7202];
C        = [0.5336   -0.2416    0.4056    0.0167    0.0115   -0.6213    0.4604    0.1259
            0.1910    0.1106    0.3421   -0.1956    0.2760    0.9206    0.4323   -0.7717
            0.2772    0.2116   -0.2094    0.4855    0.0615   -0.3314   -0.7378   -0.6914
           -0.6986   -0.0247   -0.0576    0.2415    0.3123    0.2021   -0.0033   -0.1274];
D        = zeros(4,4);
Syst{4}  = ss(A,B,C,D);
clear A B C D;

%% Turner, Kerr, and Sofrony - IJRNC 2015 - Example 17
A      = [-5.4966   -1.7013   -1.7006   -1.0871   -0.3885   -0.2560
           8.0000         0         0         0         0         0
                0    2.0000         0         0         0         0
                0         0    2.0000         0         0         0
                0         0         0    2.0000         0         0
                0         0         0         0    0.5000         0];
B       = [0.5751    0.3840    0.0158    0.6315
           0.4514    0.6831    0.0164    0.7176
           0.0439    0.0928    0.1901    0.6927
           0.0272    0.0353    0.5869    0.0841
           0.3127    0.6124    0.0576    0.4544
           0.0129    0.6085    0.3676    0.4418];
C       = [0.3533    0.7275    0.4508    0.2548    0.9084    0.0784
           0.1536    0.4784    0.7159    0.8656    0.2319    0.6408
           0.6756    0.5548    0.8928    0.2324    0.2393    0.1909
           0.6992    0.1210    0.2731    0.8049    0.0498    0.8439];
D       = zeros(4,4);
Syst{5} = ss(A,B,C,D);
clear A B C D;

%% Turner, Kerr, and Sofrony - IJRNC 2015 - Example 19
A = [-5.4966   -1.7013   -1.7006   -1.0871   -0.3885   -0.2560
      8.0000         0         0         0         0         0
           0    2.0000         0         0         0         0
           0         0    2.0000         0         0         0
           0         0         0    2.0000         0         0
           0         0         0         0    0.5000         0];
B       = [-0.1747    0.1242   -0.1618   -0.4348
            0.1069   -0.4609    0.3498    0.2012
           -0.2069   -0.2954   -0.1913    0.5389
           -0.3018    0.4648    0.6416    0.3872
           -0.3921   -0.0244    0.4655   -0.5162
           -0.8676   -0.6227    0.1863    0.2499];
C       = [-0.0874   -0.4380    0.2355   -0.1751   -0.3313   -0.0954
           -0.5379    0.0904   -0.0862    0.0422    0.0770   -0.2257
           -0.3988    0.1378   -0.1092    0.2767    0.1642   -0.2763
            0.0327   -0.4871   -0.4322    0.2185   -0.8372    0.2476];
D       = zeros(4,4);
Syst{6} = ss(A,B,C,D);
clear A B C D;

%% Turner, Kerr, and Sofrony - IJRNC 2015 - Example 23
A        = [-5.8097   -1.7258   -1.3496   -0.8882   -0.3991   -0.2843   -0.3009   -0.3524
             8.0000         0         0         0         0         0         0         0
             0    2.0000         0         0         0         0         0         0
             0         0    2.0000         0         0         0         0         0
             0         0         0    2.0000         0         0         0         0
             0         0         0         0    1.0000         0         0         0
             0         0         0         0         0    0.5000         0         0
             0         0         0         0         0         0    0.2500         0];
B        = [-0.0697    0.1982   -0.4746   -0.3289
            -0.4567    0.3788   -0.2270    0.5888
             0.3737   -0.6601    0.6961   -0.3923
            -0.7929   -0.6439   -0.0461   -0.2640
             0.0148    0.4926    0.5402    0.7878
            -0.5594   -0.0041    0.2616    0.4186
             0.0385   -0.5487   -0.5681    0.5668
            -0.1095    0.0133    0.6951    0.7202];
C        = 0.01*[0.5336   -0.2416    0.4056    0.0167    0.0115   -0.6213    0.4604    0.1259
            0.1910    0.1106    0.3421   -0.1956    0.2760    0.9206    0.4323   -0.7717
            0.2772    0.2116   -0.2094    0.4855    0.0615   -0.3314   -0.7378   -0.6914
           -0.6986   -0.0247   -0.0576    0.2415    0.3123    0.2021   -0.0033   -0.1274];
D        = zeros(4,4);
Syst{7}  = ss(A,B,C,D);
clear A B C D;

%% Drummond, Guiver, and Turner - IEEETAC 2022 - Example 2
A       = diag([-1 -10 -30 -60 -100])+0.1*ones(5,5);
B       = eye(5);
C       = 1*ones(5,5)-eye(5);
D       = zeros(5,5);
Syst{8} = ss(A,B,C,D);
clear A B C D;
