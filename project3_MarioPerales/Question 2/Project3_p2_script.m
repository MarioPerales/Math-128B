mu = .2;
X = 2;
dy = .05;
v_0 = .2;
N = 50;
theta = zeros(N,1);
max_iter = 100;

theta_2 = project3_2(theta,mu,dy,v_0,X,N,max_iter);