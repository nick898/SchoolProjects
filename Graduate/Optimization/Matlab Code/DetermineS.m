function [s] = DetermineS(c)
% This function determines the sign of each element of a vector 'c' and
% outputs the sign of each element of 'c' in a vector 's'
%
% INPUTS:
%        c - any vector
%
% OUTPUTS
%        s - a vector of signs for each element of 'c' so s is a vector 
%            containing only the values {-1,0,1}

%Initialize some variables
N = length(c);
s = zeros(1,N);

%Determine the sign of each element of c and store in s
for i = 1:N
    s(i) = sign(c(i));
end

end

