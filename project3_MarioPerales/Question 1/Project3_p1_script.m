mu = 0;
X = 2;
dy = .2;
N = 20;
v_0 = 0;
theta = ones(N,1);
max_iter = 100;

theta_1 = project3_1(theta,mu,dy,v_0,X,N,max_iter);

x = [0; cumsum(dy*tan(theta_1))]; 
y = -dy * (0:N)'; 
plot(x,y,'.-');
axis equal;
drawnow;