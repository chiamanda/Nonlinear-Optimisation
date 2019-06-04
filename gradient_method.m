% The aim of this programming assignment is to implement the gradient  descent method together with the Armijo step size rule in Matlab

function [x,grad,iter] = gradient_method(f,df,x_0)
% The function "gradient_method" receives three arguments: f and df - the function handles for the function and its gradient, and x_0 - the initial vector
% The function returns the approximate solution "x", the euclidean norm "grad" of the gradient and the number "num" of iterations.

gamma = 1; 
delta = 0.1; % flattening parameter
tol = 0.0001; % termination tolerance 
nummax = 1000; % maximum number of iterations

iter  = 0; % begin the iteration counter at 0
x = x_0; % initialize optimization vector
grad = norm(df(x)); % the Euclidean norm of the gradient of f at the point x = x_0

while (grad>=tol & iter<=nummax) % begin the gradient descent method
    d = df(x); % consider the descent direction d
    sigma = gamma *(dot(d,d) / norm(d).^2); % the algorithm begins with a sufficiently large step-size
    
    while (f(x-sigma*d) > f(x) - delta*sigma*dot(d,d)) % we want to find the first step size that satisfies the Armijo rule
        sigma = 0.5 * sigma; % the step size in each pass is halved
    end
    x = x-sigma*d; % take a step in the direction d
    
    % update while loop condition variables
    iter = iter+1; % update the iteration counter
    grad = norm(df(x)); % calculate the euclidean norm of the new x value
end
end

