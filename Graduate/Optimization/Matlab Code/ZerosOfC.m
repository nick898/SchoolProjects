function [IndexVec] = ZerosOfC(c)
% This function will find all the indices of a vector for which the
% elements of that vector are 0. It returns a vector of indices whose
% entries represent the indices of the inpute vector that are 0
%
% INPUTS:
%       c - input vector
%
% OUTPUTS: 
%       IndexVec - vector of indicies that give the indices of the input
%                  vector which are zero. If an element of IndexVec is 0
%                  that means that the corresponding element in c is
%                  non-zero

    %Initialize some variables
    N = length(c);
    IndexVec = zeros(1,N);
    
    %Determine the indices of 'c' which are zero and store in the 'IndexVec'
    for i = 1:N
        if c(i) == 0 %If the entry of 'c' is 0 then we store that indice in 
                     %'IndexVec'
            IndexVec(i) = i; %Store the indice
        else
            IndexVec(i)= 0; %'c' is non-zero so do not store the indice
        end
    end

end

