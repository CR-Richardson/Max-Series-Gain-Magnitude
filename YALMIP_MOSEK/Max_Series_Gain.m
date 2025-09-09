
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
% This script builds various linear systems and computes the maximum series
% gain (alpha) according to various criteria for which the Lurie system is
% stable when the repeated magnitude is placed in the feebdack path. 
% 
%
% Scripts
% Examples:     Contains example linear systems
%
% Functions
% LoopShift1: Maps Lurie system with ReLU nonlinearity to equivalent Lurie system with magnitude nonlinearity.
% LoopShift2: Maps Lurie system with magnitude nonlinearity to equivalent Lurie system with ReLU nonlinearity.
% Quad_Lyap:  Quadratic Criterion (Theorem 1).
% Lurie_type: Lurie-based Criterion (Corollary 2).
% SGT:        Small Gain Theorem (Equation 5.15, Boyd, 1994).
%
% Variables
% Total_Ex: Total number of examples.
% Ex_array: The set of example systems to compute alpha for.
% eps:      Loop termination accuracy - try 1e-6 for Examples 1-8 and 1e-2 for Examples 9-12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script variables
clear all; close all;
Total_Ex = 12;
Ex_array = 1:8;
eps      = 1e-6;

%% Makes example systems
Examples;

%% Nyquist gains of each example
alpha_up = [100000.0, 89.9, 0.6983, 0.0020, 0.0869, 0.8202, ...
            0.2002, 2.0221, 2.1, 2.8, 2.3, 2.8];

alpha_up = alpha_up + 0.1; % ensures alpha_low != alpha_up for each example.

%% Arrays for storing the maximum series gain (alpha) and # of decision variables.

All_ex = 1:Total_Ex; % All example systems

% Arrays for storing the maximum series gain (alpha) for each example.
alpha_Quad  = zeros(All_ex);
alpha_Lurie = zeros(All_ex);
alpha_SGT   = zeros(All_ex);

% Arrays for storing the # of decision variables for each example.
decs_Quad  = zeros(All_ex);
decs_Lurie = zeros(All_ex);
decs_SGT   = zeros(All_ex);

%% Calculate maximum series gain using various criteria

for i=Ex_array
    disp(['Example ',num2str(i),' ']);
    
    disp('Small Gain Theorem calculations...');
    tic;
    [alpha_SGT(i), data_SGT(i), decs_SGT(i)] = SGT(Syst{i}, eps, alpha_up(i));
    toc;

    disp('Quadratic Lyapunov calculations...');
    tic;
    [alpha_Quad(i), data_Quad(i), decs_Quad(i)] = Quad_Lyap(Syst{i}, eps, alpha_up(i), alpha_SGT(i));
    toc;

    disp('Lurie-type Lyapunov (H=I) calculations...');
    tic;
    [alpha_Lurie(i), data_Lurie(i), decs_Lurie(i)] = Lurie_type(Syst{i}, eps, alpha_up(i), alpha_Quad(i));
    toc;
end

%% Display max series gain

disp(' ');
disp('Max. series gain');
title_str=['        Example', '    Quadratic', '   Lurie-type', '          SGT'];
mat_data =[Ex_array' alpha_Quad(Ex_array)' alpha_Lurie(Ex_array)' alpha_SGT(Ex_array)'];
fprintf('%15s %15s %15s %15s\n',title_str);
disp(' ');
fprintf('%15d %14.4f %14.4f %14.4f\n',mat_data');

%% Display # of decision variables

disp(' ');
disp('# of decision variables');
title_str=['        Example', '    Quadratic', '   Lurie-type', '          SGT'];
mat_data =[Ex_array' decs_Quad(Ex_array)' decs_Lurie(Ex_array)' decs_SGT(Ex_array)'];
fprintf('%15s %15s %15s %15s\n',title_str);
disp(' ');
fprintf('%15d %14.0f %14.0f %14.0f\n',mat_data');
