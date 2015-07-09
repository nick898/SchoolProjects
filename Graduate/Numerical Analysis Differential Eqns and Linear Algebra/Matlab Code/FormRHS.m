function [RHS,bstar] = FormRHS(InitialV,Solution,N,D,h,j)
% This function FormRHS forms the RHS of the ODE
%   
% INPUTS: InitialV, Solution, N, D, h, j
% 
% InitialV is the initial (2 x N-1) matrix which we will use to begin all of our
% iterations to determine the actual v we use when updating the solution
% matrix for the next time step
% 
% Solution is the solution matrix. It is of size (J + 1 x N + 1) where J is
% the number of time steps and N is the number of spatial steps
% 
% N is the number of spatial steps
% 
% D is the diffusion coefficient
% 
% h is the time step 
% 
% j refers to the jth row of the Solution matrix u
% 
% OUTPUTS: RHS
% 
% RHS is the output of the function which is the RHS of the ODE. It is of
% size (2 x N-1)
% 
% 
% 

%These are the A coefficients from the Radau 2A Butcher Tableau
a = [5/12 -1/12;3/4 1/4];

%k is the spatial step
k = 1/N;

%Bu is B_{k}u which is the linear part of the ODE
Bu = diff([0 Solution(j,2:N) 0],2)*D/k/k;

%bstar is the nonlinear part of the ODE. 
bstar = phi(repmat(Solution(j,2:N),2,1)+h*a*InitialV);

%Our output, RHS, is Bu + bstar provided we replicate Bu so it is of
%appropriate dimensions to add to bstar
RHS = repmat(Bu,2,1) + bstar;

end


