function [IndexVec] = CheckC(func,xvals,yvals,gamma,c,s,parameter)
% This function finds the indices of a vector for which a certain condition
% is satisfied. The condition that needs to be met is:
%
%          |GradF(i)| > parameter AND c(i) = 0
%
% If the above condition is satisfied then we will store the indice of that
% entry. Otherwise, we set the element of the output 'IndexVec' to 0
%
% INPUTS:
%        func - function to compute the gradient of f
%        xvals - needed to evaluate the function
%        yvals - needed to evaluate the function
%        gamma - needed to evaluate the function
%        c - vector we will check for entries that are 0
%        s - needed to evaluate the function
%        parameter - this is the scalar value we will use to check if
%                    |GradF(i)| > parameter
%
% OUTPUTS
%        IndexVec - a vector storing indices of GradF for which the
%                   conditions are met

%Initialize variables
N = length(c);
IndexVec = zeros(1,N);

%Compute gradient of function
[FVal,GradF] = func(xvals,yvals,gamma,c,s,parameter);

%Find the indices of the gradient of f for which the following is true
for i = 1:N
    if (GradF(i) > parameter || GradF(i) < -parameter) && c(i) == 0
        IndexVec(i) = i; 
    else
        IndexVec(i) = 0;
    end
end


end

