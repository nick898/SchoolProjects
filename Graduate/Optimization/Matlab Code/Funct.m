function [FVal,GradF,HessF] = Funct(xvals,yvals,gamma,c,s,parameter)
% This function defines the Log-Likelihood Function by
%
%   L(gamma, c) = sum_{i = 1}^{N} [y_{i}*log(1 + exp(gamma + c'x_{i})) + ...
%                                (1 - y_{i})*log(1 + exp(c'x_{i}))]
%
% This function will also compute and return the Gradient of L and the
% Hessian Matrix of L at different values for (gamma, c)
%
% INPUTS:
%       xvals - matrix where each row is a data point with 19 features with
%               each feature representing some characteristic of a person
%       yvals - a vector of outcomes either 0 (person dies) or 1 (person
%               lives)
%       gamma - scalar value...we will look to find an optimal gamma that
%               minimizes this log-likelihood function
%       c - vector...we will look to find an optimal c that minimizes the
%           log-likelihood function
%       s - vector of signs for each element of c
%       parameter - this parameter controls how many entries of c are
%       non-zero
%
% OUTPUTS:
%       FVal - function value at a specified (gamma, c)
%       GradF - gradient of function at a specificed (gamma,c)
%       HessF - Hessian matrix of function at a specificed (gamma,c)
%

%Get size of c to determine size of hessian matrix
N = size(c,2);

%------------------------- COMPUTING THE FUNCTION VALUE --------------------


FVal = log(1 + exp(gamma + c*xvals'))*yvals + log(1 + exp(-gamma - c*xvals'))*(1 - yvals);

%Penalization
FVal = FVal + parameter*(s*c');



%------------------------- COMPUTING THE GRADIENT -------------------------
if nargout > 1
    
    GradF_g = yvals'*(exp(gamma + c*xvals')./(1 + exp(gamma + c*xvals')))' + ...
              (yvals - 1)'*(exp(-gamma - c*xvals')./(1+exp(-gamma - c*xvals')))';
     
    
    GradF_c = xvals'*((exp(gamma + c*xvals')./(1 + exp(gamma + c*xvals')))'.*yvals) + ...
              xvals'*((exp(-gamma - c*xvals')./(1 + exp(-gamma - c*xvals')))'.*(yvals - 1));
      
    %Penalization
    GradF_c = GradF_c' + parameter*s;
    
    %Define the gradient of the function as a row vector
    GradF = [GradF_c,GradF_g];
           
end



%------------------------- COMPUTING THE HESSIAN --------------------------
    
if nargout > 2
       
    HessF = zeros(N + 1, N+1);
    %Second partial w.r.t gamma (scalar)
    HessF_g = exp(gamma + c*xvals')./(1 + exp(gamma + c*xvals')).^2*yvals + ...
              exp(-gamma - c*xvals')./(1 + exp(-gamma - c*xvals')).^2*(1 - yvals);
      
   
    %Second partial w.r.t c and gamma (vector)
    HessF_cg = xvals'*(yvals'.*(exp(gamma + c*xvals')./(1 + exp(gamma + c*xvals')).^2))' + ...
               xvals'*((1 - yvals)'.*(exp(-gamma - c*xvals')./(1 + exp(-gamma - c*xvals')).^2))';
              
    %Second partial w.r.t c (matrix)
    HessF_cc = xvals'*(diag(((exp(gamma + c*xvals')./(1 + exp(gamma + c*xvals')).^2)'.*yvals)) + ...
               diag(((exp(-gamma - c*xvals')./(1 + exp(-gamma - c*xvals')).^2)'.*(1 - yvals))))*xvals; 
                  
                  
    %Define the Hessian Matrix
    HessF(end,end) = HessF_g; %Entry in last row and last column is the second partial w.r.t gamma
    HessF(end,1:end-1) = HessF_cg; %Vector in the last row up to the second to last column is the second partial w.r.t c and gamma
    HessF(1:end-1,end) = HessF_cg'; %Vector in the last column up to the second to last row is the second partial w.r.t. c and gamma
    HessF(1:end-1,1:end-1) = HessF_cc; %Matrix from element (1,1) through element (end-1,end-1) is the second partial w.r.t. c
      
end 
end