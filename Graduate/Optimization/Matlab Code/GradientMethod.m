function [FVal,c,gamma,s,GradF] = GradientMethod(func,xvals,yvals,gamma,c,s,parameter,p,IndexVec)
% This function performs a Gradient Step for the function for both gamma and
% c
%
% INPUTS:
%       func - function we want to minimize
%       xvals - needed to evaluate the function
%       yvals - needed to evaluate the function
%       gamma - want to find an optimal gamma to minimize the function
%       c - want to find an optimal c to minimize the function
%       s - vector of signs of each element of c
%       parameter - controls how many entries of c are non-zero
%       p - controls the decay rate of the Newton step size alpha
%       IndexVec - determines the entries of c which we will update with a 
%                  gradient step and the entries of c we DO NOT update with
%                  a gradient step
%
% OUTPUTS:
%       FVal - function value after taking the Gradient Step for both gamma
%              and c
%       c - updated c vector after taking Gradient Step
%       gamma - updated scalar gamma after taking Gradient Step
%       s - updated vector of signs depending on the new vector c
%       GradF - gradient of f at the new values for (gamma,c)
%

    %Initialize some variables
    N = length(c);
    GradientStep_C = zeros(1,N); %This will store the Gradient step for the vector c
    GradientStep_g = zeros(1,1); %This will store the Gradient step for gama

    %Compute initial function value, gradient, and hessian
    [FVal0,GradF,HessF] = func(xvals,yvals,gamma,c,s,parameter);

    
%------------------------------- BEGIN GRADIENT METHOD -------------------------------    
    %Set the Gradient Directions for c and gamma
    Direction_C = -GradF(1:N);
    Direction_g = -GradF(end);
    Direction = [Direction_C,Direction_g]; %This is the full gradient step as one vector

    %Compute the initial value of alpha
    alpha0 = InitialAlpha(c,Direction);

    %Do a Gradient Step ONLY for the entries of c that are either nonzero
    %or have |GradF(1:N)| > alpha
    for i = 1:N
        if IndexVec(i) ~= 0 
            GradientStep_C(i) = c(i) + alpha0*Direction_C(i)'; %Gradient step for the entry of c
        elseif c(i) ~= 0
            GradientStep_C(i) = c(i) + alpha0*Direction_C(i)'; %Gradient step for the non-zero entries of c
        else
            GradientStep_C(i) = c(i); %Keep the entry of c the same
        end
    end
    
    %Do a Gradient Step for gamma
    GradientStep_g = gamma + alpha0*Direction_g;

    %This is the full Gradient Step [d, delta]
    GradientStep = [GradientStep_C,GradientStep_g];

    %Compute the new function value after taking the Gradient Step for c and
    %gamma
    s = DetermineS(GradientStep_C); %Update s
    [FVal1,GradF,HessF] = func(xvals,yvals,GradientStep_g,GradientStep_C,s,parameter);
    alpha = alpha0; %Set alpha to be equal to the initial alpha

%------------------------------- BEGIN LINE SEARCH -------------------------------    
    while FVal1 > FVal0 
    
        alpha = p*alpha; %Decrease alpha
       
        %Do a Gradient Step ONLY for the entries of c that are either nonzero
        %or have |GradF(1:N)| > alpha
        for i = 1:N
            if IndexVec(i) ~= 0 
                GradientStep_C(i) = c(i) + alpha*Direction_C(i)'; %Gradient step for 
                                                                  %the entries of c that 
                                                                  %are zero and large in magnitude
            elseif c(i) ~= 0
                GradientStep_C(i) = c(i) + alpha*Direction_C(i)'; %Gradient step for the non-zero entries of c
            else
                GradientStep_C(i) = c(i); %Keep the entry of c the same otherwise
            end
        end
   
        %Do a Gradient Step for gamma
        GradientStep_g = gamma + alpha*Direction_g;
    
        %Compute updated function evaluation, gradient, and hessian
        s = DetermineS(GradientStep_C); %Update s
        [FVal1,GradF] = func(xvals,yvals,GradientStep_g,GradientStep_C,s,parameter);
    
    end
%------------------------------- END LINE SEARCH -------------------------------    

    %Now we have exited, the line search, we can finally update s and 
    %compute new function value and gradient
    s = DetermineS(GradientStep_C);
    [FVal1,GradF] = func(xvals,yvals,GradientStep_g,GradientStep_C,s,parameter);

%------------------------------- BEGIN NEWTONS METHOD -------------------------------    


    %Update c, s, gamma, and the function value. These will be the outputs of
    %this function along with the gradient of f
    c = GradientStep_C;
    s = DetermineS(c);
    gamma = GradientStep_g;
    FVal = FVal1;


end

