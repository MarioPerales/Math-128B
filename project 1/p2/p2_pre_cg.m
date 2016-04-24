function [x,number_of_pre_cg_iterations] = p2_pre_cg(A,b,k)

true_solution = A\b;

M = tril(A,-1) - tril(A,-2) + triu(A,1) - triu(A,2) + diag(diag(A));
M_inv = M^-1;
x = 0*b;
r = b;
p = M_inv*r;
z = p;

i = 1;

initial_error = norm(x - true_solution,inf)/norm(true_solution,inf);

while i <= k
    rdotz = r'*z;
    Ap = A*p;
    
    alpha = (rdotz)/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    z = M_inv*r;
    beta = (r'*z)/(rdotz);
    p = z + beta*p;
    
    pre_cg_error = norm(true_solution - x,inf)/norm(true_solution,inf);
    if(pre_cg_error < 10e-6*initial_error)
        break
    end
    i = i + 1;
end
number_of_pre_cg_iterations = i;
if(number_of_pre_cg_iterations > 1000);
    number_of_pre_cg_iterations = 0;
end