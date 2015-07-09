function [xvals_N] = Normalize(xvals)
% This function will normalize the 1st, 9th, and 10th columns of the xvals
% matrix. Since the values in these three columns are so different than the
% other values in the matrix then normalizing will allow us to do
% computations so that function values do not blow up
%
% INPUTS:   
%       xvals - matrix of data points
%
% OUTPUTS:
%       xvals_N - matrix where the 1st, 9th, and 10th columns have been
%                 normalized by dividing the infinity norm of the column
%

%Define the normalized xvals matrix
xvals_N = xvals;

%Normalize the 1st, 9th, and 10th columns of xvals_N by dividing by the
%infinity norm of the respective column
xvals_N(:,1) = xvals_N(:,1)/norm(xvals_N(:,1),inf);
xvals_N(:,9) = xvals_N(:,9)/norm(xvals_N(:,9),inf);
xvals_N(:,10) = xvals_N(:,10)/norm(xvals_N(:,10),inf);



end

