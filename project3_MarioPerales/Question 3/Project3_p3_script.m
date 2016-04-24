mu = .2;
X = 2;
dy = .02;
v_0 = .2;
N = 50;
max_iter = 100;
theta = zeros(N,1);
N_c = 4; % Change this for varying results...

theta_3 = project3_3(theta,mu,dy,v_0,X,N,max_iter,N_c);

