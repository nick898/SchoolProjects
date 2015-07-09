%README:
%
%       ONLY NEED TO HAVE THE 'data' MATRIX FROM THE FILE 'read_icu_data_file' IN
%       YOUR MATLAB WORKSPACE. YOU CAN OBTAIN THE 'read_icu_data_file' FROM Icon
%
%       ALTERNATIVELY, YOU COULD ALSO JUST LOAD THE 'DataWorkspace' THAT I
%       HAVE INCLUDED IN THE ZIP FILE
%
%       We can change all of the various parameters below
%
%

%Necessary Data (DO NOT CHANGE)
xvals = data(:,3:end);
xvals_N = Normalize(xvals);
yvals = data(:,2);

%DEFINE THE VARIOUS PARAMETERS (FEEL FREE TO CHANGE)
tol = 1e-10; %Tolerance threshold
c = zeros(1,19); %Initial vector for c
gamma = 0; %Initial scalar value for gamma
parameter = 0.05; %This parameter controls the number of non-zero entries in c
s = DetermineS(c); %Determine the signs of each element of c and store this in the vector s
p = 0.5; %Define the decay rate of the step size in the line search

%Counter to keep track of the number of iterations
Counter = 0;

%Compute the gradient of the lagrangian function to begin the iterations
%Grad_Lagrangian = GradL(@Funct,parameter,gamma,c,xvals_N,yvals,s);

FVals = zeros(1,50);
GradientCell = cell(1,50);
HessianCell = cell(1,50);
GradLagrangeCell = cell(1,50);
GradFuncNorm = zeros(1,50);
GradLagrangeNorm = zeros(1,50);


tic
Grad_Lagrangian = GradL(@Funct,parameter,gamma,c,xvals_N,yvals,s);

while norm(Grad_Lagrangian) > tol    
    
    %Newton Method
    [FVal,c,gamma,s,GradF] = NewtonMethod(@Funct,xvals_N,yvals,gamma,c,s,parameter,p);
        
    %Check to see if the partial derivative w.r.t. c is small enough
    [IndexVec] = CheckC(@Funct,xvals_N,yvals,gamma,c,s,parameter);
        
    %Gradient Method
    [FVal,c,gamma,s,GradF] = GradientMethod(@Funct,xvals_N,yvals,gamma,c,s,parameter,p,IndexVec);
  
    %Store the function values, gradient vectors, and hessian matrices for
    %each iteration in a cell array
    FVals(Counter + 1) = FVal;
    GradientCell{Counter + 1} = GradF;
    GradFuncNorm(Counter + 1) = norm(GradF);
    
    %Update the gradient of the Lagrangian
    Grad_Lagrangian = GradL(@Funct,parameter,gamma,c,xvals_N,yvals,s);
    GradLagrangeCell{Counter + 1} = Grad_Lagrangian;
    GradLagrangeNorm(Counter + 1) = norm(Grad_Lagrangian);
    
    %Keep a counter on the number of iterations
    Counter = Counter + 1;

end
toc

[Grad_Lagrangian,Lambda_m,Lambda_p] = GradL(@Funct,parameter,gamma,c,xvals,yvals,s);
















