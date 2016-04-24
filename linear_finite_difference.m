function linear_finite_difference(p,q,r,a_given,b_given,alpha,beta,N,exact_solution)

if(nargin < 9)
    exact_solution = @(x) 0;
end

x_vector = zeros(N+2,1);

w = zeros(N+2,1);

h = (b_given-a_given)/(N+1);
x = a_given + h;

a = zeros(N+1,1); %a
b = zeros(N+1,1); %b
c = zeros(N+1,1); %c
d = zeros(N+1,1); %d
l = zeros(N+1,1); %l
u = zeros(N+1,1); %u
z = zeros(N+1,1); %z

a(2) = 2 + h^2 * q(x);
b(2) = -1 + (h/2)*p(x);
d(2) = -h^2*r(x) + (1 + (h/2)*p(x))*alpha;

for i = 3:N
    x = a_given + (i-1)*h;
    a(i) = 2 + h^2*q(x);
    b(i) = -1 + (h/2)*p(x);
    c(i) = -1 - (h/2)*p(x);
    d(i) = -h^2*r(x);
end
x = b_given-h;
a(N+1) = 2 + h^2*q(x);
c(N+1) = -1 - (h/2)*p(x);
d(N+1) = -h^2*r(x) + (1 - (h/2)*p(x))*beta;

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

k = N;
w(1) = alpha;
w(N+2) = beta;
w(N+1) = z(N+1);
while k >= 2;
    w(k) = z(k) - u(k)*w(k+1);
    k = k - 1;
end

fprintf('\n\nw_1 ~ y(x)\n\n');
fprintf('      x              w_1         |y(x) - w_1|\n');
fprintf('----------------------------------------------\n');
for i = 1:N+2;
    x = a_given + (i-1)*h;
    x_vector(i) = x;
    fprintf('%12.9f    %12.9f    %12.9f \n',x,w(i),abs(exact_solution(x) - w(i)));
end

