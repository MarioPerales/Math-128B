function [x, pre_cg_error] = pre_cg(A,b,k)

true_solution = A\b;
pre_cg_error = zeros(1,k);

M = tril(A,-1) - tril(A,-2) + triu(A,1) - triu(A,2) + diag(diag(A));
M_inv = M^-1;
x = 0*b;
r = b;
p = M_inv*r;
z = p;

i = 1;
while i <= k
    rdotz = r'*z;
    Ap = A*p;
    
    alpha = (rdotz)/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    z = M_inv*r;
    beta = (r'*z)/(rdotz);
    p = z + beta*p;
    
    pre_cg_error(i) = norm(true_solution - x,inf)/norm(true_solution,inf);
    i = i + 1;
end
    