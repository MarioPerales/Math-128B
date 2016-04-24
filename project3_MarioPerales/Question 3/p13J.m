function J = p13J(theta,mu,dy,v_0,X,N)
J = zeros(N,N);
delta = 1e-6;
e_i = zeros(N,1);

for i = 1:N
    e_i(i) = 1;
    J(:,i) = (p13F(theta + e_i*delta,mu,dy,v_0,X,N) - p13F(theta - e_i*delta,mu,dy,v_0,X,N))/(2*delta);
    e_i(i) = 0;
end