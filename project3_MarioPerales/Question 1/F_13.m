function F = F_13(theta,mu,dy,X,v_0,N)

F = zeros(N,1);
v = zeros(N,1);
w = zeros(N,1);
g = 32.17;

for j = 1:N
    v(j) = sqrt(v_0^2 + 2*g*j*dy - 2*mu*dy*sum(1./cos(theta(1:j))));
end

for k = 1:N
    w(k) = -dy*v(j)*sum(1./(v(1:N).^3.*cos(theta(1:N))));
end

for i = 1:N-1
    F(i) = sin(theta(i+1))/v(i+1)*(1 - mu*w(i+1)) - sin(theta(i))/v(i)*(1 - mu*w(i));
end
F(N) = dy*sum(tan(theta)) - X;