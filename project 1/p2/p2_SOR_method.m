function [x, number_of_SOR_iterations] = p2_SOR_method(A,b,omega,k)
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
initial_error = norm(x - true_solution,inf)/norm(true_solution,inf);

while i <= k
    x = T*x+c;
    SOR_error = norm(x - true_solution,inf)/norm(true_solution,inf);
    if(SOR_error < 10e-6*initial_error)
        break
    end
    i = i + 1;
end
number_of_SOR_iterations = i;
if(number_of_SOR_iterations > 1000);
    number_of_SOR_iterations = 0;
end