function [x, cg_error] =my_cg(A,b,k)

true_solution = A\b;

cg_error = zeros(1,k);
x = 0*b;
r = b;
v = b;

i = 1;

while i <= k
    rdotr = r'*r;
    Av = A*v;
    
    t = (rdotr)/(v'*Av);
    x = x + t*v;
    r = r - t*Av;
    s = (r'*r)/(rdotr);
    v = r + s*v;
    
    cg_error(i) = norm(true_solution - x,inf)/norm(true_solution,inf);
    i = i + 1;
end