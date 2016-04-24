function [c,h] = piecewise_linear_rayleigh_ritz(p,q,f,x)

n = length(x) - 2;

h = zeros(n+1,1);
Q1 = zeros(n+1,1);
Q2 = zeros(n+1,1);
Q3 = zeros(n+1,1);
Q4 = zeros(n+2,1);
Q5 = zeros(n+1,1);
Q6 = zeros(n+1,1);
alpha = zeros(n+1,1);
beta = zeros(n+1,1);
a = zeros(n+1,1);
b = zeros(n+1,1);
c = zeros(n+1,1);
chi = zeros(n+1,1);
z = zeros(n+1,1);

for i=1:n+1
    h(i) = x(i+1) - x(i);
end

for i = 2:n
    Q1(i) = (1/h(i))^2*integral(@(y) (x(i+1)-y).*(y - x(i)).*q(y),x(i),x(i+1));
    Q2(i) = (1/h(i-1))^2*integral(@(y) (y-x(i-1)).^2.*q(y),x(i-1),x(i));
    Q3(i) = (1/h(i))^2*integral(@(y) (x(i+1) - y).^2.*q(y),x(i),x(i+1));
    Q4(i) = (1/h(i-1))^2*integral(p,x(i-1),x(i));
    Q5(i) = 1/h(i-1)*integral(@(y) (y-x(i-1)).*f(y),x(i-1),x(i));
    Q6(i) = 1/h(i)*integral(@(y) (x(i+1)-y).*f(y),x(i),x(i+1));
end

Q1(n+1) = (1/h(n+1))^2*integral(@(y) (x(n+2)-y).*(y - x(n+1)).*q(y),x(n+1),x(n+2));
Q2(n+1) =(1/h(n))^2*integral(@(y) (y-x(n)).^2.*q(y),x(n),x(n+1));
Q3(n+1) = (1/h(n+1))^2*integral(@(y) (x(n+2) - y).^2.*q(y),x(n+1),x(i+2));
Q4(n+1) = (1/h(n))^2*integral(p,x(n),x(n+1));
Q4(n+2) =(1/h(n+1))^2*integral(p,x(n+1),x(n+2));
Q5(n+1) = 1/h(n)*integral(@(y) (y-x(n)).*f(y),x(n),x(n+1));
Q6(n+1) = 1/h(n+1)*integral(@(y) (x(n+2)-y).*f(y),x(n+1),x(n+2));

for i=2:n
    alpha(i) = Q4(i) + Q4(i+1) + Q2(i) + Q3(i);
    beta(i) = Q1(i) - Q4(i+1);
    b(i) = Q5(i) + Q6(i);
end

alpha(n+1) = Q4(n+1) + Q4(n+2) + Q2(n+1) + Q3(n+1);
b(n+1) = Q5(n+1) + Q6(n+1);

a(2) = alpha(2);
chi(2) = beta(2)/alpha(2);
z(2) = b(2)/a(2);

for i = 3:n
    a(i) = alpha(i) - beta(i-1)*chi(i-1);
    chi(i) = beta(i)/a(i);
    z(i) = (b(i) - beta(i-1) * z(i-1))/a(i);
end

a(n+1) = alpha(n+1) - beta(n)*chi(n);
z(n+1) = (b(n+1) - beta(n)*z(n))/a(n+1);

c(n+1) = z(n+1);

fprintf('c%d = %12.9f\n',n,c(n+1));

k = n;


while k >= 2;
    c(k) = z(k) - chi(k)*c(k+1);
    fprintf('c%d = %12.9f\n',k-1,c(k));
    k = k-1;
end

c = c(2:n+1);
end

