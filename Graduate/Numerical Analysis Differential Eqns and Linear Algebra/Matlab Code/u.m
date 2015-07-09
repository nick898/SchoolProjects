function [Solution] = u(t,x)
% This function computes our solution matrix
% 
% INPUTS t,x
% 
% t is the time domain which should always be [0, 20] at steps of h
% 
% x is the spatial domain which should always be [0,1] at steps of k = 1/N
% where N is the number of spatial steps
% 
% OUTPUTS: Solution
% 
% Solution will be a matrix which will be our numerical approximation. It 
% is of size (J + 1 x N + 1) where J is the number of time steps and N is 
% the number of spatial steps
%
%

% Spatial Discretization
N = length(x) - 1; %Number of spatial steps
% k = 1/N; Spatial step. I don't use this, but I like to see it just for
% reference

% Time Discretization...I don't use J or h in this program this, but I like 
% to see it just for reference
% J = length(t) - 1; Number of time steps
% h = 20/J; Time Step

%Initialize our Solution as a matrix
Solution = zeros(length(t),length(x));
Solution(:,1) = 0; %Set Boundary condition
Solution(:,end) = 0; %Set Boundary Condition

%%%%%%%%%%%%%%%  VARIOUS INITIAL CONDITIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial Condition 1:
%Set the entire first column of our Solution. So for x = 0, we know the
%solution for all time. This is the test initial value function that
%Stewart asked us to use in the assignment
for i = 1:N-1
    Solution(1,i+1) = sin(3*pi*i/N);
end

%Initial Condition 2:
%Random Initial Value Function
% Solution(1,2:100) = rand(1,N-1) - rand(1,N-1); 

%Initial Condition 3:
%Random Initial Value Function
% Solution(1,2:100) = randn(1,N-1);

%Initial Condition 4:
% Other Known Initial Value
% for i = 1:N-1
%     Solution(1,i+1) = sin(2*pi*i/N) - sin(5*pi*i/N);
% end

%Initial Condition 5:
% %Initial Condition all 1's
% Solution(1,2:N) = 1;



end

