# Max-Series-Gain-Magnitude
This code acompanies the paper [Analysis of Lurie Systems with Magnitude Nonlinearities and Connections to Neural Network Stability Analysis]() where the experimental setup is detailed in the *numerical examples* section. The code computes the maximum series gain for which global asymptotic stability is verified using two of the criteria presented in the paper. Furthermore, the number of decision variables used by the criteria is also returned to compare the complexity. The criteria are tested on a number of example Lurie systems assumed to have repeated ReLU nonlinearities, where the relevant loopshift is used to allow the criteria below to be applied. Implementations of the other criteria can be found in the repository [Max-Series-Gain](https://github.com/CR-Richardson/Max-Series-Gain).

### Authors:
* Carl R Richardson (cr2g16@soton.ac.uk)
* Matthew C Turner (m.c.turner@soton.ac.uk)

## Prerequisites
All the code is written in MATLAB. The LMI's are solved using the *Robust Control Toolbox* which must be installed as an add-on.

## Overview
The repository is organised as follows:
- `Max_Series_Gain.m` The master script. It loops through each example, computing the maximum series gain (and # of decision variables) according to each criterion,  and displays the results.
- `Examples.m` Defines the (A,B,C,D) matrices of the example systems.
- `LoopShift1.m` Performs the loopshift from Lurie system with ReLU nonlinearity to Lurie system with magnitude nonlinearity.
- `LoopShift2.m` Performs the loopshift from Lurie system with magnitude nonlinearity to Lurie system with ReLU nonlinearity.
- `Quad_Lyap.m` Implementation of the Quadratic Criterion - See Theorem 1.
- `Lurie_type.m` Implementation of the Lurie-based Criterion - See Corollary 2.
- `Aizerman.m` Computes the Nyquist gain (based on the Aizerman Conjecture) for each example.

## Getting Started
Run `Max_Series_Gain.m` to repeat the experiments in the paper or select a subset of the examples by defining them in the *Ex_array* variable.  
