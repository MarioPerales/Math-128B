function [x, jacobi_error] = jacobi_method(A,b,k)
jacobi_error = zeros(1,k);

matrix_size = size(A);
N = matrix_size(1);
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
x = zeros(N,1);
i = 1;
T = D^-1*(L+U);
c = D^-1*b;

true_solution = A\b;

while i <= k
    x = T*x+c;
    jacobi_error(i) = norm(x - true_solution,inf)/norm(true_solution,inf);
    
    i = i + 1;
end