function [Solution] = USolver(InitialV,Solution,A,B,N,D,h)
% This is the main program which we will use to compute the numerical
% approximation to the Reaction-Diffusion Equation. Running this program
% directly will call all the other necessary functions needed to compute
% the approximation to the RD equation
% 
% INPUTS: InitialV, Solution, A, B, N, D, h
% 
% InitialV is the initial (2 x N-1) matrix which we will use to begin all of our
% iterations to determine the actual V we use when updating the solution
% matrix for the next time step
%
% Solution is the Solution matrix which our numerical approximation will
% be stored in. It is of size (J + 1 x N + 1) where J is the number of time 
% steps and N is the number of spatial steps
%
% A (2x2 matrix) are the blocks along the main diagonal of the larger
% matrix block tridiaginal matrix
%
% B (2x2 matrix)are the blocks along both the sub-diagonal and
% super-diagonal of the larger block tridiagonal matrix
% 
% N is the number of spatial steps
% 
% D is the Diffusion Coefficient
% 
% h is the time step
%
% OUTPUTS: Solution
%
% Solution is the Solution matrix. It is of size (J + 1 x N + 1) where J is
% the number of time steps and N is the number of spatial steps
%

 for j = 1:1000 %length(Solution) - 1
    v = IterativeSolver(InitialV,Solution,A,B,N,D,h,j);


Solution(j+1,2:N) = Solution(j,2:N) + h*(3/4*v(1,:) + 1/4*v(2,:));
end

end

