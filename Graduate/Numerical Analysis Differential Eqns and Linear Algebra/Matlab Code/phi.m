function [b] = phi(u)
%This function defines phi, which is the nonlinear part of the ODE
% 
% INPUTS: u
%
% u is our solution matrix
% 
% OUTPUTS: b
% 
% b is phi(u)
%
%


% Test Phi that Stewart asked us to use in the assignment
b = u.*(1 - u.^2);

% Setting Phi = 0...this is a useful way to test our solver to see if what
% we get is the solution to the Heat Eqn
% b = 0;

end

