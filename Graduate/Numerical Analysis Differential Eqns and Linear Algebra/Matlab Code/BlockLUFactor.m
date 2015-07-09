function [v] = BlockLUFactor(A,B,RHS)
%This function computes the LU Factorization of a Block Tridiagonal Matrix. 
%
%INPUTS: A (2x2 matrix), B (2x2 matrix), and RHS (2 x N-1 matrix)
%
%The A (2x2 matrix) are the blocks along the main diagonal of the larger
%matrix block tridiaginal matrix
%
%The B (2x2 matrix)are the blocks along both the sub-diagonal and
%super-diagonal of the larger block tridiagonal matrix
%
%The RHS is the RHS of the ODE. It is a (2 x N-1) matrix. (In my comments,
%when I refer to m, just remember m = N-1
%
%
%OUTPUTS: v
%
%v is the solution after doing the LU decomposition and doing the forward/backward
%substitutions
%
%

%Set m (m = N-1)
m = length(RHS);

%Reshape RHS so it's a single column
RHS = reshape(RHS,2*m, 1);


%Define L,M,U, V
L = cell(m,1);
M = cell(m,1);
U = cell(m,1);
V = cell(m,1);

%Set L{1},U{1},V{1}, and M{1}
[L{1},U{1}] = lu(A); 
V{1} = L{1}\B;
M{1} = B/U{1};

%Loop to find L{k}, U{k}, V{k}, and M{k} for k = 2,...,m
for k = 2:m
    [L{k},U{k}] = lu(A - M{k-1}*V{k-1});
    V{k} = L{k}\B;
    M{k} = B/U{k};
end

%Forward Substitution
y = zeros(2*m,1); %Define y

y(1:2,:) = L{1}\RHS(1:2,:); %Define the first y

%Loop for Forward Substitution
for k = 1:m-1
    y(2*k+1:2*k+2,:) = L{k+1}\(RHS(2*k+1:2*k+2,:) - M{k}*y(2*k-1:2*k,:));
end

%Back Substitution
x = zeros(2*m,1); %Define x

x(2*m-1:2*m,:) = U{m}\y(2*m-1:2*m,:); %Know the last x

%Loop for Backward Substitution
for k = m-1:-1:1
    x(2*k-1:2*k,:) = U{k}\(y(2*k-1:2*k,:) - V{k}*x(2*k+1:2*k+2,:));
end

%Turn x into a (2xm matrix) and define this to be our output: v (Remember m = N-1)
v = reshape(x,2,m);

end

