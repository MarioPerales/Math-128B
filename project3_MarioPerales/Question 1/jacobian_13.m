function J = jacobian_13(theta,mu,dy,v_0,N)
J = zeros(N,N);
v = zeros(N,1);
g = 32.17;

for j = 1:N
    v(j) = sqrt(v_0^2 + 2*g*j*dy - 2*mu*dy*sum(1./cos(theta(1:j))));
end

for i = 1:N-1
    J(i,i) = -cos(theta(i))/v(i);
    J(i,i+1) = cos(theta(i+1))/v(i+1);
end

for k = 1:N
    J(N,k) = dy*(sec(theta(k)))^2;
end
