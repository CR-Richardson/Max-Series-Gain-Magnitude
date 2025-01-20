
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
% This script builds various linear systems and computes the maximum series
% gain (alpha) according to various criteria for which the Lurie system is
% stable when the repeated magnitude is placed in the feebdack path. 
% 
% Note: Only Quadratic Criterion can handle D~=0.
%
% Scripts
% Examples:     Contains example linear systems
%
% Functions
% LoopShift1: Maps Lurie system with ReLU nonlinearity to equivalent Lurie system with magnitude nonlinearity.
% LoopShift2: Maps Lurie system with magnitude nonlinearity to equivalent Lurie system with ReLU nonlinearity.
% Quad_Lyap:  Quadratic Criterion (Theorem 1).
% Lurie_type: Lurie-based Criterion (Corollary 2).

%
% Variables
% Total_Ex: Total number of examples
% Ex_array: The set of example systems to compute alpha for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script variables
clear all; close all;
Total_Ex   = 8;
Ex_array   = 1:8;

%% Makes example systems
Examples;

%% Calculate maximum series gain using various criteria

All_ex = 1:Total_Ex; % All example systems

% Arrays for storing the maximum series gain (alpha) for each example.
Alpha_Quad = zeros(All_ex);
Alpha_Lurie = zeros(All_ex);

% Arrays for storing the # of decision variables for each example.
decs_Quad = zeros(All_ex);
decs_Lurie = zeros(All_ex);

for i=Ex_arrayli
    disp(['Example ',num2str(i),' ']);
    
    disp('Quadratic Lyapunov calculations...'); 
    [Alpha_Quad(i), data1(i), decs_Quad(i)] = Quad_Lyap(Syst{i});

    disp('Lurie-type Lyapunov (H=I) calculations...'); 
    [Alpha_Lurie(i), data2(i), decs_Lurie(i)] = Lurie_type(Syst{i});
end

%% Display max series gain

disp(' ');
disp('Max. series gain');
title_str=['        Example', '    Quadratic', '   Lurie-type'];
mat_data =[Ex_array' Alpha_Quad(Ex_array)' Alpha_Lurie(Ex_array)'];
fprintf('%15s %15s %15s\n',title_str);
disp(' ');
fprintf('%15d %14.4f %14.4f\n',mat_data');

%% Display # of decision variables

disp(' ');
disp('# of decision variables');
title_str=['        Example', '    Quadratic', '   Lurie-type'];
mat_data =[Ex_array' decs_Quad(Ex_array)' decs_Lurie(Ex_array)'];
fprintf('%15s %15s %15s\n',title_str);
disp(' ');
fprintf('%15d %14.0f %14.0f\n',mat_data');
