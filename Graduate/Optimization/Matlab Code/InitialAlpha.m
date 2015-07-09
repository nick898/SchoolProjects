function [alpha0] = InitialAlpha(c,d)
% This function computes an initial step size alpha0 which satisfies a
% certain condition stated below:
%
%   alpha0 = min[ 1, min_{c>0,d<0} (-c_{i}/d_{i}), min_{c<0, d>0} (+c_{i}/d_{i}) ]
%
% This will allow us to conduct a line search limited by the condition that
% the sign of each c_{i} does not change
%
% INPUTS:
%       c - vector c
%       d - Newton direction vector for c
%
% OUTPUTS: 
%       alpha0 - initial step size alpha to begin line search

%Initialize some variables
n = size(c,2);
a1 = 1;
a2 = 1;

%Compute min_{c>0,d<0} (-c_{i}/d_{i})
for i = 1:n
    if c(i) > 0 && d(i) < 0
        Temp = -c(i)/d(i);
        if a1 > Temp
            a1 = Temp;
        else
            a1 = a1;
        end
    end
end

%Compute min_{c<0, d>0} (+c_{i}/d_{i})
for i = 1:n
    if c(i) < 0 && d(i) > 0
        Temp = -c(i)/d(i);
        if a2 > Temp
            a2 = Temp;
        else
            a2 = a2;
        end
    end
end

%Set the initial alpha as being the minimum between a1 and a2
alpha0 = min(a1,a2);

end