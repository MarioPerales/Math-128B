function [x, number_of_pre_cg_ichol_iterations, pre_cg_ichol_error] = p2_pre_cg_ichol(A,b,k)

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
initial_error = norm(x - true_solution,inf)/norm(true_solution,inf);

while i <= k
    rdotz = r'*z;
    Ap = A*p;
    
    alpha = (rdotz)/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    z = R\(R'\r);
    beta = (r'*z)/(rdotz);
    p = z + beta*p;
    
    pre_cg_ichol_error = norm(true_solution - x,inf)/norm(true_solution,inf);
    if(pre_cg_ichol_error < 10e-6*initial_error)
        break
    end
    i = i + 1;
end
number_of_pre_cg_ichol_iterations = i;
if(number_of_pre_cg_ichol_iterations > 1000);
    number_of_pre_cg_ichol_iterations = 0;
end    