function [x, pre_cg_ichol_error] = pre_cg_ichol(A,b,k)

true_solution = A\b;
pre_cg_ichol_error = zeros(1,k);

R = ichol(A);
% M = R'*R; %% need to fix
% M_inv = M^-1;
x = 0*b;
r = b;
p = R\(R'\r);
z = p;

i = 1;
while i <= k
    rdotz = r'*z;
    Ap = A*p;
    
    alpha = (rdotz)/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    z = R\(R'\r);
    beta = (r'*z)/(rdotz);
    p = z + beta*p;
    
    pre_cg_ichol_error(i) = norm(true_solution - x,inf)/norm(true_solution,inf);
    i = i + 1;
end
    