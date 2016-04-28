%------------------------------%
%          Constants           %
%------------------------------%

n1 = 5;
n2 = 10;
n3 = 20;
n4 = 40;
n5 = 80;
max_iter = 1000;
z = @(x,y) 3.*y.*(1-y);
sol_n1 = zeros(n1+1,n1+1);
sol_n2 = zeros(n2+1,n2+1);
sol_n3 = zeros(n3+1,n3+1);
sol_n4 = zeros(n4+1,n4+1);
sol_n5 = zeros(n5+1,n5+1);

%------------------------------%
%         First Grid           %
%------------------------------%

h = 1/n1;
k = 1;
[x,y] = ndgrid(0:h:1, 0:h:1);
flag = 0;
A = speye((n1+1)^2, (n1+1)^2);

fprintf('\n\n\n n = %d\n\n\n', n1)
fprintf('k         norm(z_k - z_k-1)_inf\n')
fprintf('-------------------------------\n')

while k <= max_iter
    n = n1;
    b = z(x(:), y(:));
    for i = 2:n
        for j = 2:n
            
            z_x = (sol_n1(i+1,j) - sol_n1(i-1,j))/(2*h);
            z_y = (sol_n1(i,j+1) - sol_n1(i,j-1))/(2*h);
            c1 = 1+z_x^2;
            c2 = 2*z_x*z_y;
            c3 = 1+z_y^2;
            
            t = i + (j-1)*(n+1);
            A(t,t) = -2*c1/h^2 + -2*c3/h^2; %#ok<*SPRIX>
            A(t,t+1) = c3/h^2;
            A(t,t-1) = c3/h^2;
            A(t,t + (n+1)) = c1/h^2;
            A(t,t - (n+1)) = c1/h^2;
            A(t,(t+1) + (n+1)) = -c2/(4*h^2);
            A(t,(t+1) - (n+1)) = c2/(4*h^2);
            A(t,(t-1) + (n+1)) = c2/(4*h^2);
            A(t,(t-1) - (n+1)) = -c2/(4*h^2);
            
            b(t) = 0;
        end
    end
    new_sol_n1 = reshape(A\b,size(x));
    new_sol_vec_n1 = A\b;
    old_sol_vec_n1 = reshape(sol_n1,[(n+1)^2,1]);
    
    fprintf('%d         %12.8f\n',k+1,norm(new_sol_vec_n1 - old_sol_vec_n1,inf))
    if(norm(new_sol_vec_n1 - old_sol_vec_n1,inf) < 1e-8)
        sol_n1 = new_sol_n1;
        flag = 1;
        break;
    end
    sol_n1 = new_sol_n1;
    k = k + 1;
end

if(flag == 0)
    fprintf('\n\n\n No Convergence')
end


%------------------------------%
%        Second Grid           %
%------------------------------%

h = 1/n2;
k = 1;
[x,y] = ndgrid(0:h:1, 0:h:1);
flag = 0;
A = speye((n2+1)^2, (n2+1)^2);

fprintf('\n\n\n n = %d\n\n\n', n2)
fprintf('k         norm(z_k - z_k-1)_inf\n')
fprintf('-------------------------------\n')

while k <= max_iter
    n = n2;
    b = z(x(:), y(:));
    for i = 2:n
        for j = 2:n
            
            z_x = (sol_n2(i+1,j) - sol_n2(i-1,j))/(2*h);
            z_y = (sol_n2(i,j+1) - sol_n2(i,j-1))/(2*h);
            c1 = 1+z_x^2;
            c2 = 2*z_x*z_y;
            c3 = 1+z_y^2;
            
            t = i + (j-1)*(n+1);
            A(t,t) = -2*c1/h^2 + -2*c3/h^2; %#ok<*SPRIX>
            A(t,t+1) = c3/h^2;
            A(t,t-1) = c3/h^2;
            A(t,t + (n+1)) = c1/h^2;
            A(t,t - (n+1)) = c1/h^2;
            A(t,(t+1) + (n+1)) = -c2/(4*h^2);
            A(t,(t+1) - (n+1)) = c2/(4*h^2);
            A(t,(t-1) + (n+1)) = c2/(4*h^2);
            A(t,(t-1) - (n+1)) = -c2/(4*h^2);
            
            b(t) = 0;
        end
    end
    new_sol_n2 = reshape(A\b,size(x));
    new_sol_vec_n2 = A\b;
    old_sol_vec_n2 = reshape(sol_n2,[(n+1)^2,1]);
    
    fprintf('%d         %12.8f\n',k+1,norm(new_sol_vec_n2 - old_sol_vec_n2,inf))
    if(norm(new_sol_vec_n2 - old_sol_vec_n2,inf) < 1e-8)
        sol_n2 = new_sol_n2;
        flag = 1;
        break;
    end
    sol_n2 = new_sol_n2;
    k = k + 1;
end

if(flag == 0)
    fprintf('\n\n\n No Convergence')
end

%------------------------------%
%         Third Grid           %
%------------------------------%

h = 1/n3;
k = 1;
[x,y] = ndgrid(0:h:1, 0:h:1);
flag = 0;
A = speye((n3+1)^2, (n3+1)^2);

fprintf('\n\n\n n = %d\n\n\n', n3)
fprintf('k         norm(z_k - z_k-1)_inf\n')
fprintf('-------------------------------\n')

while k <= max_iter
    n = n3;
    b = z(x(:), y(:));
    for i = 2:n
        for j = 2:n
            
            z_x = (sol_n3(i+1,j) - sol_n3(i-1,j))/(2*h);
            z_y = (sol_n3(i,j+1) - sol_n3(i,j-1))/(2*h);
            c1 = 1+z_x^2;
            c2 = 2*z_x*z_y;
            c3 = 1+z_y^2;
            
            t = i + (j-1)*(n+1);
            A(t,t) = -2*c1/h^2 + -2*c3/h^2; %#ok<*SPRIX>
            A(t,t+1) = c3/h^2;
            A(t,t-1) = c3/h^2;
            A(t,t + (n+1)) = c1/h^2;
            A(t,t - (n+1)) = c1/h^2;
            A(t,(t+1) + (n+1)) = -c2/(4*h^2);
            A(t,(t+1) - (n+1)) = c2/(4*h^2);
            A(t,(t-1) + (n+1)) = c2/(4*h^2);
            A(t,(t-1) - (n+1)) = -c2/(4*h^2);
            
            b(t) = 0;
        end
    end
    new_sol_n3 = reshape(A\b,size(x));
    new_sol_vec_n3 = A\b;
    old_sol_vec_n3 = reshape(sol_n3,[(n+1)^2,1]);
    
    fprintf('%d         %12.8f\n',k+1,norm(new_sol_vec_n3 - old_sol_vec_n3,inf))
    if(norm(new_sol_vec_n3 - old_sol_vec_n3,inf) < 1e-8)
        sol_n3 = new_sol_n3;
        flag = 1;
        break;
    end
    sol_n3 = new_sol_n3;
    k = k + 1;
end

if(flag == 0)
    fprintf('\n\n\n No Convergence')
end

%------------------------------%
%        Fourth Grid           %
%------------------------------%

h = 1/n4;
k = 1;
[x,y] = ndgrid(0:h:1, 0:h:1);
flag = 0;
A = speye((n4+1)^2, (n4+1)^2);

fprintf('\n\n\n n = %d\n\n\n', n4)
fprintf('k         norm(z_k - z_k-1)_inf\n')
fprintf('-------------------------------\n')

while k <= max_iter
    n = n4;
    b = z(x(:), y(:));
    for i = 2:n
        for j = 2:n
            
            z_x = (sol_n4(i+1,j) - sol_n4(i-1,j))/(2*h);
            z_y = (sol_n4(i,j+1) - sol_n4(i,j-1))/(2*h);
            c1 = 1+z_x^2;
            c2 = 2*z_x*z_y;
            c3 = 1+z_y^2;
            
            t = i + (j-1)*(n+1);
            A(t,t) = -2*c1/h^2 + -2*c3/h^2; %#ok<*SPRIX>
            A(t,t+1) = c3/h^2;
            A(t,t-1) = c3/h^2;
            A(t,t + (n+1)) = c1/h^2;
            A(t,t - (n+1)) = c1/h^2;
            A(t,(t+1) + (n+1)) = -c2/(4*h^2);
            A(t,(t+1) - (n+1)) = c2/(4*h^2);
            A(t,(t-1) + (n+1)) = c2/(4*h^2);
            A(t,(t-1) - (n+1)) = -c2/(4*h^2);
            
            b(t) = 0;
        end
    end
    new_sol_n4 = reshape(A\b,size(x));
    new_sol_vec_n4 = A\b;
    old_sol_vec_n4 = reshape(sol_n4,[(n+1)^2,1]);
    
    fprintf('%d         %12.8f\n',k+1,norm(new_sol_vec_n4 - old_sol_vec_n4,inf))
    if(norm(new_sol_vec_n4 - old_sol_vec_n4,inf) < 1e-8)
        sol_n4 = new_sol_n4;
        flag = 1;
        break;
    end
    sol_n4 = new_sol_n4;
    k = k + 1;
end

if(flag == 0)
    fprintf('\n\n\n No Convergence')
end

%------------------------------%
%         Fifth Grid           %
%------------------------------%

h = 1/n5;
k = 1;
[x,y] = ndgrid(0:h:1, 0:h:1);
flag = 0;
A = speye((n5+1)^2, (n5+1)^2);

fprintf('\n\n\n n = %d\n\n\n', n5)
fprintf('k         norm(z_k - z_k-1)_inf\n')
fprintf('-------------------------------\n')

while k <= max_iter
    n = n5;
    b = z(x(:), y(:));
    for i = 2:n
        for j = 2:n
            
            z_x = (sol_n5(i+1,j) - sol_n5(i-1,j))/(2*h);
            z_y = (sol_n5(i,j+1) - sol_n5(i,j-1))/(2*h);
            c1 = 1+z_x^2;
            c2 = 2*z_x*z_y;
            c3 = 1+z_y^2;
            
            t = i + (j-1)*(n+1);
            A(t,t) = -2*c1/h^2 + -2*c3/h^2; %#ok<*SPRIX>
            A(t,t+1) = c3/h^2;
            A(t,t-1) = c3/h^2;
            A(t,t + (n+1)) = c1/h^2;
            A(t,t - (n+1)) = c1/h^2;
            A(t,(t+1) + (n+1)) = -c2/(4*h^2);
            A(t,(t+1) - (n+1)) = c2/(4*h^2);
            A(t,(t-1) + (n+1)) = c2/(4*h^2);
            A(t,(t-1) - (n+1)) = -c2/(4*h^2);
            
            b(t) = 0;
        end
    end
    new_sol_n5 = reshape(A\b,size(x));
    new_sol_vec_n5 = A\b;
    old_sol_vec_n5 = reshape(sol_n5,[(n+1)^2,1]);
    
    fprintf('%d         %12.8f\n',k+1,norm(new_sol_vec_n5 - old_sol_vec_n5,inf))
    if(norm(new_sol_vec_n5 - old_sol_vec_n5,inf) < 1e-8)
        sol_n5 = new_sol_n5;
        flag = 1;
        break;
    end
    sol_n5 = new_sol_n5;
    k = k + 1;
end

if(flag == 0)
    fprintf('\n\n\n No Convergence')
end

real_sol = sol_n5;

% We create a set of mesh points for each solution we find to compare it 
% to the real solution. We collect the correct mesh points in the following
% lines of code.

%------------------------------%
%      First mesh points       %
%------------------------------%

real_sol_n1_mesh_points = zeros((n1+1)^2,1);
steps_n1 = n5/n1;
u = 1;
col_iter = 1;
k = 1;

while k <= (n1+1)^2
    row_iter = mod((k-1)*steps_n1,n5+steps_n1) + 1;
    if(row_iter == n5 + 1)
        real_sol_n1_mesh_points(k) = real_sol(row_iter,col_iter);
        col_iter = u*steps_n1 + 1;
        u = u + 1;
        k = k + 1;
        continue
    end
    real_sol_n1_mesh_points(k) = real_sol(row_iter,col_iter);
    k = k + 1;
end

%------------------------------%
%      Second mesh points      %
%------------------------------%


real_sol_n2_mesh_points = zeros((n2+1)^2,1);
steps_n2 = n5/n2;
u = 1;
col_iter = 1;
k = 1;

while k <= (n2+1)^2
    row_iter = mod((k-1)*steps_n2,n5+steps_n2) + 1;
    if(row_iter == n5 + 1)
        real_sol_n2_mesh_points(k) = real_sol(row_iter,col_iter);
        col_iter = u*steps_n2 + 1;
        u = u + 1;
        k = k + 1;
        continue
    end
    real_sol_n2_mesh_points(k) = real_sol(row_iter,col_iter);
    k = k + 1;
end

%------------------------------%
%       Third mesh points      %
%------------------------------%

real_sol_n3_mesh_points = zeros((n3+1)^2,1);
steps_n3 = n5/n3;
u = 1;
col_iter = 1;
k = 1;

while k <= (n3+1)^2
    row_iter = mod((k-1)*steps_n3,n5+steps_n3) + 1;
    if(row_iter == n5 + 1)
        real_sol_n3_mesh_points(k) = real_sol(row_iter,col_iter);
        col_iter = u*steps_n3 + 1;
        u = u + 1;
        k = k + 1;
        continue
    end
    real_sol_n3_mesh_points(k) = real_sol(row_iter,col_iter);
    k = k + 1;
end

%------------------------------%
%      Fourth mesh points      %
%------------------------------%

real_sol_n4_mesh_points = zeros((n4+1)^2,1);
steps_n4 = n5/n4;
u = 1;
col_iter = 1;
k = 1;

while k <= (n4+1)^2
    row_iter = mod((k-1)*steps_n4,n5+steps_n4) + 1;
    if(row_iter == n5 + 1)
        real_sol_n4_mesh_points(k) = real_sol(row_iter,col_iter);
        col_iter = u*steps_n4 + 1;
        u = u + 1;
        k = k + 1;
        continue
    end
    real_sol_n4_mesh_points(k) = real_sol(row_iter,col_iter);
    k = k + 1;
end

% Compare real solution mesh points to each corresponding solution.

n1_error = norm(real_sol_n1_mesh_points - new_sol_vec_n1,inf);
n2_error = norm(real_sol_n2_mesh_points - new_sol_vec_n2,inf);
n3_error = norm(real_sol_n3_mesh_points - new_sol_vec_n3,inf);
n4_error = norm(real_sol_n4_mesh_points - new_sol_vec_n4,inf);

% Plot

y = [n1_error;n2_error;n3_error;n4_error];
x = [1/n1;1/n2;1/n3;1/n4];

loglog(x,y)
xlabel('h');
ylabel('Max norm error');
title('Problem 2');

% Convergence Rate

coeff = polyfit(x,y,1);
rate = exp(coeff(1));

fprintf('\n\n\n The convergence rate is%12.8f\n', rate)
