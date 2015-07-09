%This Script can be used to easily change certain parameters and quickly
%enter them in Matlab as we will need them to run the code and we are
%constantly changing them. 
% 
% This script will automatically generate various constants that are needed 
% as well as the matrices A and B that form the Tridiagonal Block Matrix
% 
% In addition to this it will initialize the space and time domains as well
% as the random vector InitialV which we will use to begin the Iterative
% Solver. It also initializes the Solution matrix through the function
% u(t,x). This Solution matrix that is generated can be entered directly as
% is into USolver and USolver will compute the numerical approximation to 
% the Reaction-Diffusion Equation.
% 
% So all you need to do is call this script initially and then go to
% USolver and copy and paste the declaration statement in USolver into the
% Matlab Command Window.
% 
% The programs that are part of this project are listed below with a short
% summary of what they do:
% 
% 1. BlockLUFactor - LU Factorization of a Block Tridiagonal Matrix
% 2. IterativeSolver - Iterative Solver to find V1 and V2 used to update
%                      the (n+1)'st time step
% 3. FormRHS - Forms the RHS of the PDE to be used in IterativeSolver
% 4. USolver - the main program which computes the numerical approximation.
%              This calls all the necessary functions to compute the
%              approximation.
% 5. phi - this forms the nonlinear part of the PDE
% 6. u - this forms the Solution matrix initially
% 


%This sets the constants D,N,k, and h
D = 10^(-2); %Set Diffusion Coefficient
N = 100; %Set N, the number of spatial steps
%Other Value of N
%N = 200;
k = 1/N; %Set k, the spatial step
h = 0.01; %Set h, the time step
%Other value of h
% h = .005;


%This sets the matrices A and B which are part of the tridiagonal Block
%Matrix
A = eye(2) + 2*h*D*[5/12 -1/12;3/4 1/4]/k/k; %This is the block matrix along the main diagonal
B = -h*D*[5/12 -1/12;3/4 1/4]/k/k; %This is the block matrix along the sub-diagonals and super-diagonals

%Define space and time intervals
x = 0:k:1; %Spatial domain
t = 0:h:20; %Time domain

%Initialize Random Vector V to run iterations
% InitialV = randn(2,N-1); %This is the initial random vector v we will use to run the iteration

%Initialize V as a matrix of zeroes if you'd like
InitialV = zeros(2,N-1);

%Initialize Solution Matrix
Solution = u(t,x);

