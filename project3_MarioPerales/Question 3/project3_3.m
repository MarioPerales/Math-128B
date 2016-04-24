function theta = project3_3(theta,mu,dy,v_0,X,N,max_iter,N_c)

h = 1/N_c;
b = -h*p13F(theta,mu,dy,v_0,X,N);

for i = 1:N_c
    A = p13J(theta,mu,dy,v_0,X,N);
    k1 = A\b;
    A = p13J(theta + 1/2*k1,mu,dy,v_0,X,N);
    k2 = A\b;
    A = p13J(theta + 1/2*k2,mu,dy,v_0,X,N);
    k3 = A\b;
    A = p13J(theta + k3,mu,dy,v_0,X,N);
    k4 = A\b;
    
    theta = theta + (k1 + 2*k2 + 2*k3 + k4)/6;
end

theta = project3_2(theta,mu,dy,v_0,X,N,max_iter);