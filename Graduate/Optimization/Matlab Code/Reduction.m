function [Grad,Hess] = Reduction(s,Grad,Hess)
% This function will reduce the size of the gradient of f and the size of
% the hessian matrix of f depending on whether or not an entry of the
% vector 's' is zero
%
% INPUTS: 
%       s - vector containing the signs of each element of the vector c
%       Grad - gradient vector of f
%       Hess - hessian matrix of f
%
% OUTPUTS:
%       Grad - gradient vector reduced in size depending on the entries of
%              's' that are 0
%       Hess - hessian matrix reduced in size depending on the entries of
%              's' that are 0
%
%

%Initialize a variable
N = length(s);

%If an element of the vector 's' is 0, then we ignore the corresponding
%gradient entry and ignore the corresponding row and column of the hessian
%matrix
for i = N:-1:1
    if s(i) == 0
        Grad(i) = [];
        Hess(i,:) = [];
        Hess(:,i) = [];
    end
end



end

