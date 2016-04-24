function [x, SOR_error] = SOR_method(A,b,omega,k)

SOR_error = zeros(1,k);
matrix_size = size(A);
N = matrix_size(1);
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
x = zeros(N,1);
i = 1;

T = (D-omega*L)^-1*((1-omega)*D + omega*U);
c = omega*(D-omega*L)^-1*b;

true_solution = A\b;

while i <= k
    x = T*x+c;
    SOR_error(i) = norm(x - true_solution,inf)/norm(true_solution,inf);
    
    i = i + 1;
end