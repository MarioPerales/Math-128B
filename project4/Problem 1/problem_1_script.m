%------------------------------%
%          Constants           %
%------------------------------%

n = 40;
h = 1/n;
max_iter = 1000;
k = 1;
[x,y] = ndgrid(0:h:1, 0:h:1);
z = @(x,y) 3.*y.*(1-y);
sol = zeros(n+1,n+1);
flag = 0;
A = speye((n+1)^2, (n+1)^2);

fprintf('k         norm(z_k - z_k-1)_inf\n')
fprintf('-------------------------------\n')

%------------------------------%
%           Main Loop          %
%------------------------------%

while k <= max_iter
    b = z(x(:), y(:));
    for i = 2:n
        for j = 2:n
            
            z_x = (sol(i+1,j) - sol(i-1,j))/(2*h);
            z_y = (sol(i,j+1) - sol(i,j-1))/(2*h);
            c1 = 1+z_x^2;
            c2 = 2*z_x*z_y;
            c3 = 1+z_y^2;
            
            t = i + (j-1)*(n+1);
            
            A(t,t) = -2*c1/h^2 + -2*c3/h^2; %#ok<*SPRIX>
            A(t,t+1) = c3/h^2; % z(i+1,j)
            A(t,t-1) = c3/h^2; % z(i-1,j)
            A(t,t + (n+1)) = c1/h^2; % z(i,j+1)
            A(t,t - (n+1)) = c1/h^2; % z(i,j-1)
            A(t,(t+1) + (n+1)) = -c2/(4*h^2); % z(i+1,j+1)
            A(t,(t+1) - (n+1)) = c2/(4*h^2); % z(i+1,j-1)
            A(t,(t-1) + (n+1)) = c2/(4*h^2); % z(i-1,j+1)
            A(t,(t-1) - (n+1)) = -c2/(4*h^2); % z(i-1,j-1)
            
            b(t) = 0;
        end
    end
    new_sol_vec = A\b;
    new_sol = reshape(new_sol_vec,size(x));
    old_sol_vec = reshape(sol,[(n+1)^2,1]);
    
    fprintf('%d         %12.8f\n',k+1,norm(new_sol_vec - old_sol_vec,inf))
    if(norm(new_sol_vec - old_sol_vec,inf) < 1e-8)
        flag = 1;
        break;
    end
    sol = new_sol;
    k = k + 1;
end

% Plot if it converged

if(flag == 1)
    surf(x', y', new_sol'), axis equal, shading interp
    cameramenu, lighting phong, camlight right
else
    fprintf('No Convergence')
end
    