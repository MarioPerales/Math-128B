function nonlinear_finite_difference(f,dfdy,dfdyp,a_given,b_given,alpha,beta,N,M,tol,exact_solution)

if(nargin < 11)
    exact_solution = @(x) 0;
end

x_vector = zeros(N+2,1);

h = (b_given-a_given)/(N+1);

w = zeros(N+2,1);
w(1) = alpha;
w(N+2) = beta;

a = zeros(N+1,1); %a
b = zeros(N+1,1); %b
c = zeros(N+1,1); %c
d = zeros(N+1,1); %d
l = zeros(N+1,1); %l
u = zeros(N+1,1); %u
v = zeros(N+1,1);
z = zeros(N+1,1); %z

for i = 2:N+1
    w(i) = alpha + (i - 1)*(beta - alpha)/(b_given - a_given)*h;
end

j = 1;

while j <= M
    x = a_given + h;
    t = (w(3) - alpha)/(2*h);
    a(2) = 2 + h^2 * dfdy(x,w(2),t);
    b(2) = -1 + (h/2)*dfdyp(x,w(2),t);
    d(2) = -(2*w(2) - w(3) - alpha + h^2*f(x,w(2),t)); 

    for i = 3:N
        x = a_given + (i-1)*h;
        t = (w(i + 1) - w(i-1))/(2*h);
        a(i) = 2 + h^2*dfdy(x,w(i),t);
        b(i) = -1 + (h/2)*dfdyp(x,w(i),t);
        c(i) = -1 - (h/2)*dfdyp(x,w(i),t);
        d(i) = -(2*w(i) - w(i+1) - w(i-1) + h^2*f(x,w(i),t));
    end
    x = b_given-h;
    t = (beta - w(N))/(2*h);
    a(N+1) = 2 + h^2*dfdy(x,w(N+1),t);
    c(N+1) = -1 - (h/2)*dfdyp(x,w(N+1),t);
    d(N+1) = -(2*w(N+1) - w(N) - beta + h^2*f(x,w(N+1),t));

    l(2) = a(2);
    u(2) = b(2)/a(2);
    z(2) = d(2)/l(2);

    for i = 3:N
        l(i) = a(i) - c(i)*u(i-1);
        u(i) = b(i)/l(i);
        z(i) = (d(i) - c(i)*z(i-1))/l(i);
    end
    l(N+1) = a(N+1) - c(N+1)*u(N);
    z(N+1) = (d(N+1) - c(N+1)*z(N))/l(N+1);

    s = N;
    v(N+1) = z(N+1);
    w(N+1) = w(N+1) + v(N+1);
    while s >= 2;
        v(s) = z(s) - u(s)*v(s+1);
        w(s) = w(s) + v(s);
        s = s - 1;
    end
    if(norm(v(2:N+1),inf) < tol)
        fprintf('\n\nw_1 ~ y(x)\n\n');
        fprintf('      x              w_1         |y(x) - w_1|\n');
        fprintf('----------------------------------------------\n');
        for i = 1:N+2;
            x = a_given + (i-1)*h;
            x_vector(i) = x;
            fprintf('%12.9f    %12.9f    %12.9f \n',x,w(i),abs(exact_solution(x) - w(i)));
        end
        break
    end
    j = j + 1;
end
