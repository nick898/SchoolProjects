function [NewV] = IterativeSolver(InitialV,Solution,A,B,N,D,h,j)
% This function is an Iterative Solver for the Block Tridiagonal Matrix
% 
% INPUTS: InitialV, Solution, A, B, N, D, h, j
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
% j refers to the jth row of the Solution matrix
% 
% OUTPUTS: NewV
%
% NewV is the (2 x N-1) matrix which we will use to update the (j+1)'st
% time step
%
%

%Form RHS initially with InitialV
RHS = FormRHS(InitialV,Solution,N,D,h,j);

%Set OldV = InitialV so we can compute the error in the first iteration
OldV = InitialV;

%Initialize error as 1
Error = 1;

%Set Epsilon
Epsilon = 10^(-10);

%While loop to run the iteration
while Error > Epsilon
    NewV = BlockLUFactor(A,B,RHS); %NewV comes out of BlockLUFactor
    Error = norm(NewV - OldV); %Update error
    OldV = NewV; %Store the just computed NewV as OldV
    RHS = FormRHS(OldV,Solution,N,D,h,j); %Form RHS again using the NewV we
                                          %just computed and enter it back 
                                          %into the loop

end

