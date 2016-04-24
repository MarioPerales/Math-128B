function [x, number_of_cg_iterations] = p2_my_cg(A,b,k)
true_solution = A\b;

x = 0*b;
r = b;
v = b;

i = 1;

initial_error = norm(x - true_solution,inf)/norm(true_solution,inf);

while i <= k
    rdotr = r'*r;
    Av = A*v;
    
    t = (rdotr)/(v'*Av);
    x = x + t*v;
    r = r - t*Av;
    s = (r'*r)/(rdotr);
    v = r + s*v;
    
    cg_error = norm(true_solution - x,inf)/norm(true_solution,inf);
    if(cg_error < 10e-6*initial_error)
        break
    end
    i = i + 1;
end
number_of_cg_iterations = i;
if(number_of_cg_iterations > 1000);
    number_of_cg_iterations = 0;
end