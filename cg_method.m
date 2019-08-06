% In this programming assignment I will implement the CG method for solving the linear system Ax=b, A:Rn×n, b in Rn
% This solution uses the notes provided in the Nonlinear Optimisation exercise on 28.05.2019

function [x,grad,num] = cg_method(A,b,x_0)
% The function "cg_method" inputs a positive definite symmetric matrix
% Abx_0 and returns the solution "x" to the linear system Ax = b

tol = 0.0001; % some fixed tolerance/precision
nummax = 1000; % a maximum number of iterations "nummax"
  
% Delete the curly bracket "{" to test the example below
%{ 
    % Example: inputs to test the function
    n = 1000;
    Q = orth(randn(n));
    D = diag(abs(randn(n,1)) + 0.5);
    A = Q*D*Q'; % symmetric positive definite matrix A
    b = rand(n,1)
    x_0 = rand(n,1)
    xsol = linsolve(A,b)
%}

% f = @(x) 0.5 * (x' *A *x)- (b' * x) % to solve the linear system we will minimise the function f (here it is not necessary to compute it) 
grad = @(x) A * x - b; % the residual computed by the gradient of the function f
   
d = -grad(x_0); % consider the initial descent direction d as the antigradient of the function f at the initial point
x =  x_0; % set x to the initial point x_0
num = 0; % begin the iteration counter

while (norm(grad(x)) > tol && num <= nummax) % Stopping criteria for the program are when the Euclidean norm of the residual is less than the tolerance or when the maximum number of iterations are reached.
    
    sigma = d' * d ./ (d' * A *d); % use the exact step size at each step
    xk = x + sigma * d; % take a step
    
    beta = (dot(-grad(x),d))/(dot(A*d,d)) % find beta 
    dk = -grad(xk) + beta * d; % the direction computed using the CG-method
    
    % update condition parameters for the solution approximation and the iteration counter, and store dk for the next step
    x = xk;    
    num = num + 1;
    d = dk;

end    