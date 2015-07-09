function [Grad_Lagrangian,Lambda_m,Lambda_p] = GradL(func,parameter,gamma,c,xvals,yvals,s)
% This function computes the Gradient of the Lagrangian of the function we
% are trying to minimize using contrained optimization
%
% INPUTS:
%       func - function we want to minimize
%       parameter - scalar valued controls the number of elements of c that
%                   are non-zero
%       gamma - scalar value...we are trying to determine an optimal gamma
%               that minimizes the function
%       c - vector...we are trying to determine an optimal c that minimizes
%           the function
%       xvals - needed to evaluate the function
%       yvals - needed to evaluate the function
%       s - vector of signs of each element in c
%
% OUTPUTS: 
%       Grad_Lagrangian - gradient of the Lagrangian function associated
%       with our constrained optimization problem for our function
%		Lambda_m - Set of Lagrange Multipliers for one constraint function
%		Lambda_p - Set of Lagrange Multipliers for the other constraint function
%



%Compute function value and gradient of function at (gamma, c)
[FVal,GradF] = func(xvals,yvals,gamma,c,s,parameter);

%Define the first component of the gradient of the lagrangian as the
%partial of f w.r.t gamma
GradL1 = GradF(end);

%Initialize some variables/storage
n = size(GradF,2);
Lambda_m = zeros(1,n); %Lambda values for one of the inequality constraint
Lambda_p  = zeros(1,n); %Lambda values for the second inequality constraint
GradL2 = zeros(1,n); %Second component of the gradient of the lagrangian
GradL3 = zeros(1,n); %Third component of the gradient of the lagrangian

%Compute the Lambdas and the second and third componenets of the gradient
%of the lagrangian
for i = 1:n
    Lambda_m(i) = 0.5*(parameter - GradF(i));
    Lambda_p(i) = 0.5*(parameter + GradF(i));
    GradL2(i) = GradF(i) + Lambda_m(i) - Lambda_p(i);
    GradL3(i) = parameter - Lambda_m(i) - Lambda_p(i);
end
    

%Define the Gradient of the Lagrangian with the three different
%componenents
Grad_Lagrangian = [GradL1, GradL2, GradL3];
end

