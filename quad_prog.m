% This programming assignment deals with minimising the function 
% f(x) = 0.5 <Qx,x> + <q,x> subject to the constraint Ax = b where Q is an
% nxn symmetric positive deinite matrix in the reals, q an arbitrary vector
% in Rn, A a full rank matrix of size mxn in the reals, where m <= n and b
% in Rn.
% This solution uses Chappter 6.1.1 of the lecture notes 

function x = solve_quad_prog(Q,q,A,b)
% The function "solve_quad_prog" inputs a positive definite symmetric matrix
% "Q", an arbitray vector ""q, a full rank matrix "A", a vector "b" and returns the solution "x" to the linear system Ax = b
  
% Delete the curly bracket "{" to test the example below
%{
    % Sample inputs to test the function
    n = 300; %  choose size n

    P = orth(randn(n)); % random orthogornal matrix
    D = diag(abs(randn(n,1))); % diagonal matrix
    Q = P*D*P'; % positive definite matrix Q
    Q=(Q+Q')/2;% symmetrise Q
    q = rand(n,1);

    m = 200;
    A = rand(m,m)*rand(m,n); %

    while (rank(A)<m) % check that A has full rank
        A = rand(m,m)*rand(m,n);
    end;
    b = rand(m,1);
%}
    ineqA = []  % set the inequality contraints to zero
    ineqb = []  % set the inequality contraints to zero
    x = quadprog(Q,q,ineqA, ineqb, A, b) % minimizes 1/2*xT*Q*x + qT*x subject to the equality constraint  A*x = b
    
% sources used: https://ch.mathworks.com/help/optim/ug/quadprog.html
    
   