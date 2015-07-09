function [FVal,c,gamma,s,GradF] = NewtonMethod(func,xvals,yvals,gamma,c,s,parameter,p)
% This function performs a Newton Step for the function for both gamma and
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
%
% OUTPUTS:
%       FVal - function value after taking the Newton Step for both gamma
%              and c
%       c - updated c vector after taking Newton Step
%       gamma - updated scalar gamma after taking Newton Step
%       s - updated vector of signs depending on the new vector c
%       GradF - gradient of f at the new values for (gamma,c)
%

    %Initialize variables
    N = length(c);
    NewtonStep_C = zeros(1,N); %This will store the Newton step for the vector c
    NewtonStep_g = zeros(1,1); %This will store the Newton step for gamma

    %Compute initial function value, gradient, and hessian
    [FVal0,GradF,HessF] = func(xvals,yvals,gamma,c,s,parameter);

    %Find the entries of c that are 0
    [IndexVec] = ZerosOfC(c);
    
    %If the entry of c is 0, then we ignore the corresponding entry in the
    %gradient and the corresponding row/column in the hessian matrix
    [GradF,HessF] = Reduction(s,GradF,HessF);
    
%------------------------------- BEGIN NEWTONS METHOD -------------------------------    
    
    %Compute the Newton direction for c and gamma
    Direction = -HessF\GradF'; 
    N2 = length(Direction);
    
    %If the length of the direction vector is less than 2 then we will return the
    %following results and exit the function
    if N2 < 2
        FVal = FVal0;
        c = c;
        s = DetermineS(c);
        alpha = InitialAlpha(c,Direction);
        return
    end
    
    %If the length of the direction vector greater than 2 then we will do a
    %Newton step for some entries of c. However, since in each iteration we
    %are ignoring some entries of the gradient and ignoring rows/columns of
    %the hessian this means that at each iteration the direction vector
    %could be a different size. So we need to account for this change in
    %size of the direction vector
    Direction_C = Direction(1:N2-1); %Newton direction for c
    Direction_C(N2:N) = 0;
    Direction_g = Direction(end); %Newton direction for gamma
    
    %This rearranges the Newton Directions for the elements of
    %c based on whether the elements of c are zero or non-zero.
    for i = 1:N
        if IndexVec(i) ~= 0
            Direction_C(i+1:N) = Direction_C(i:N-1);
            Direction_C(i) = 0;
        end
    end
    
    %Finally, this is our full Newton Direction
    Direction = [Direction_C',Direction_g];
    
    %Compute the initial value of alpha
    alpha0 = InitialAlpha(c,Direction);
    
    %Do a Newton Step ONLY for the entries of c that are non-zero
    for i = 1:N
        if IndexVec(i) == 0 
            NewtonStep_C(i) = c(i) + alpha0*Direction_C(i)'; %Newton step for the entry of c
        else
            NewtonStep_C(i) = c(i); %Keep the entry of c the same
        end
    end

    %Do a Newton Step for gamma
    NewtonStep_g = gamma + alpha0*Direction_g;

    %This is the full Newton Step [d, delta]
    NewtonStep = [NewtonStep_C,NewtonStep_g];

    %Compute the new function value after taking the Newton step for c and
    %gamma
    s = DetermineS(NewtonStep_C); %Update s 
    [FVal1,GradF,HessF] = func(xvals,yvals,NewtonStep_g,NewtonStep_C,s,parameter);
    alpha = alpha0; %Set alpha to be equal to the initial alpha
    

%------------------------------- BEGIN LINE SEARCH -------------------------------    
    while FVal1 > FVal0 
    
        alpha = p*alpha; %Decrease alpha
       
        %Do a Newton Step ONLY for the entries of c that are 0
        for i = 1:N
            if IndexVec(i) == 0
                NewtonStep_C(i) = c(i) + alpha*Direction_C(i)';
            else
                NewtonStep_C(i) = c(i);
            end
        end
   
        %Do a Newton Step for gamma
        NewtonStep_g = gamma + alpha*Direction_g;
    
        %Compute updated function evaluation, gradient, and hessian
        s = DetermineS(NewtonStep_C); %Update s
        [FVal1] = func(xvals,yvals,NewtonStep_g,NewtonStep_C,s,parameter);
    
    end
%------------------------------- END LINE SEARCH -------------------------------    
    

    %Now we have exited, the line search, we can finally update s and 
    %compute new function value and gradient
    s = DetermineS(NewtonStep_C);
    [FVal1,GradF] = func(xvals,yvals,NewtonStep_g,NewtonStep_C,s,parameter);
    
%------------------------------- END NEWTONS METHOD ------------------------------- 


    %Update c, s, gamma, and the function value. These will be the outputs of
    %this function along with the gradient of f
    c = NewtonStep_C;
    s = DetermineS(c);
    gamma = NewtonStep_g;
    FVal = FVal1;

end

