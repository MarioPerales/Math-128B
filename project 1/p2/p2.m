tic;

% You could easily forloop this but...%

k = 1000;
n1 = 5;
N1 = n1^2;
n2 = 10;
N2 = n2^2;
n3 = 20;
N3 = n3^2;
n4 = 50;
N4 = n4^2;
n5 = 100;
N5 = n5^2;

omega1 = 2/(1+sqrt(8)/(n1+1));
omega2 = 2/(1+sqrt(8)/(n2+1));
omega3 = 2/(1+sqrt(8)/(n3+1));
omega4 = 2/(1+sqrt(8)/(n4+1));
omega5 = 2/(1+sqrt(8)/(n5+1));

A1 = delsq(numgrid('S', n1+2)) * (n1+1)^2;
A2 = delsq(numgrid('S', n2+2)) * (n2+1)^2;
A3 = delsq(numgrid('S', n3+2)) * (n3+1)^2;
A4 = delsq(numgrid('S', n4+2)) * (n4+1)^2;
A5 = delsq(numgrid('S', n5+2)) * (n5+1)^2;
b1 = ones(n1^2,1);
b2 = ones(n2^2,1);
b3 = ones(n3^2,1);
b4 = ones(n4^2,1);
b5 = ones(n5^2,1);

[~, number_of_jacobi_iterations1] = p2_jacobi_method(A1,b1,k);
[~, number_of_gauss_iterations1] = p2_gauss_method(A1,b1,k);
[~, number_of_SOR_iterations1] = p2_SOR_method(A1,b1,omega1,k);
[~, number_of_cg_iterations1] = p2_my_cg(A1,b1,k);
[~, number_of_pre_cg_iterations1] = p2_pre_cg(A1,b1,k);
[~, number_of_pre_cg_ichol_iterations1] = p2_pre_cg_ichol(A1,b1,k);

[~, number_of_jacobi_iterations2] = p2_jacobi_method(A2,b2,k);
[~, number_of_gauss_iterations2] = p2_gauss_method(A2,b2,k);
[~, number_of_SOR_iterations2] = p2_SOR_method(A2,b2,omega2,k);
[~, number_of_cg_iterations2] = p2_my_cg(A2,b2,k);
[~, number_of_pre_cg_iterations2] = p2_pre_cg(A2,b2,k);
[~, number_of_pre_cg_ichol_iterations2] = p2_pre_cg_ichol(A2,b2,k);

[~, number_of_jacobi_iterations3] = p2_jacobi_method(A3,b3,k);
[~, number_of_gauss_iterations3] = p2_gauss_method(A3,b3,k);
[~, number_of_SOR_iterations3] = p2_SOR_method(A3,b3,omega3,k);
[~, number_of_cg_iterations3] = p2_my_cg(A3,b3,k);
[~, number_of_pre_cg_iterations3] = p2_pre_cg(A3,b3,k);
[~, number_of_pre_cg_ichol_iterations3] = p2_pre_cg_ichol(A3,b3,k);

[~, number_of_jacobi_iterations4] = p2_jacobi_method(A4,b4,k);
[~, number_of_gauss_iterations4] = p2_gauss_method(A4,b4,k);
[~, number_of_SOR_iterations4] = p2_SOR_method(A4,b4,omega4,k);
[~, number_of_cg_iterations4] = p2_my_cg(A4,b4,k);
[~, number_of_pre_cg_iterations4] = p2_pre_cg(A4,b4,k);
[~, number_of_pre_cg_ichol_iterations4] = p2_pre_cg_ichol(A4,b4,k);

[~, number_of_jacobi_iterations5] = p2_jacobi_method(A5,b5,k);
[~, number_of_gauss_iterations5] = p2_gauss_method(A5,b5,k);
[~, number_of_SOR_iterations5] = p2_SOR_method(A5,b5,omega5,k);
[~, number_of_cg_iterations5] = p2_my_cg(A5,b5,k);
[~, number_of_pre_cg_iterations5] = p2_pre_cg(A5,b5,k);
[~, number_of_pre_cg_ichol_iterations5] = p2_pre_cg_ichol(A5,b5,k);

N = [N1,N2,N3,N4,N5];
number_of_jacobi_iterations = [number_of_jacobi_iterations1,number_of_jacobi_iterations2,number_of_jacobi_iterations3,number_of_jacobi_iterations4,number_of_jacobi_iterations5];
number_of_gauss_iterations = [number_of_gauss_iterations1,number_of_gauss_iterations2,number_of_gauss_iterations3,number_of_gauss_iterations4,number_of_gauss_iterations5];
number_of_SOR_iterations = [number_of_SOR_iterations1,number_of_SOR_iterations2,number_of_SOR_iterations3,number_of_SOR_iterations4,number_of_SOR_iterations5];
number_of_cg_iterations = [number_of_cg_iterations1,number_of_cg_iterations2,number_of_cg_iterations3,number_of_cg_iterations4,number_of_cg_iterations5];
number_of_pre_cg_iterations = [number_of_pre_cg_iterations1,number_of_pre_cg_iterations2,number_of_pre_cg_iterations3,number_of_pre_cg_iterations4,number_of_pre_cg_iterations5];
number_of_pre_cg_ichol_iterations = [number_of_pre_cg_ichol_iterations1,number_of_pre_cg_ichol_iterations2,number_of_pre_cg_ichol_iterations3,number_of_pre_cg_ichol_iterations4,number_of_pre_cg_ichol_iterations5];

loglog(N, number_of_jacobi_iterations,'r');
hold on
loglog(N, number_of_gauss_iterations,'g')
hold on
loglog(N, number_of_SOR_iterations,'b');
hold on
loglog(N, number_of_cg_iterations,'k');
hold on
loglog(N, number_of_pre_cg_iterations,'m');
hold on
loglog(N, number_of_pre_cg_ichol_iterations,'c');

% If you don't want line plots and just plot data points, then uncomment this
% section and comment out the previous loglog plot.

% %% loglog(N1, number_of_jacobi_iterations1,'r.');
% hold on
% loglog(N1, number_of_gauss_iterations1,'g.');
% hold on
% loglog(N1, number_of_SOR_iterations1,'b.');
% hold on
% loglog(N1, number_of_cg_iterations1,'k.');
% hold on
% loglog(N1, number_of_pre_cg_iterations1,'m.');
% hold on
% loglog(N1, number_of_pre_cg_ichol_iterations1,'c.');
% hold on
% loglog(N2, number_of_jacobi_iterations2,'r.');
% hold on
% loglog(N2, number_of_gauss_iterations2,'g.');
% hold on
% loglog(N2, number_of_SOR_iterations2,'b.');
% hold on
% loglog(N2, number_of_cg_iterations2,'k.');
% hold on
% loglog(N2, number_of_pre_cg_iterations2,'m.');
% hold on
% loglog(N2, number_of_pre_cg_ichol_iterations2,'c.');
% hold on
% loglog(N3, number_of_jacobi_iterations3,'r.');
% hold on
% loglog(N3, number_of_gauss_iterations3,'g.');
% hold on
% loglog(N3, number_of_SOR_iterations3,'b.');
% hold on
% loglog(N3, number_of_cg_iterations3,'k.');
% hold on
% loglog(N3, number_of_pre_cg_iterations3,'m.');
% hold on
% loglog(N3, number_of_pre_cg_ichol_iterations3,'c.');
% hold on
% loglog(N4, number_of_jacobi_iterations4,'r.');
% hold on
% loglog(N4, number_of_gauss_iterations4,'g.');
% hold on
% loglog(N4, number_of_SOR_iterations4,'b.');
% hold on
% loglog(N4, number_of_cg_iterations4,'k.');
% hold on
% loglog(N4, number_of_pre_cg_iterations4,'m.');
% hold on
% loglog(N4, number_of_pre_cg_ichol_iterations4,'c.');
% hold on
% loglog(N5, number_of_jacobi_iterations5,'r.');
% hold on
% loglog(N5, number_of_gauss_iterations5,'g.');
% hold on
% loglog(N5, number_of_SOR_iterations5,'b.');
% hold on
% loglog(N5, number_of_cg_iterations5,'k.');
% hold on
% loglog(N5, number_of_pre_cg_iterations5,'m.');
% hold on
% loglog(N5, number_of_pre_cg_ichol_iterations5,'c.');

legend('Jacobi','Gauss-Seidel','SOR','CG', 'Pre-CG', 'Pre-CG (with ichol function)','Location','southeast');
xlabel('Size of System');
ylabel('Number of iterations');
time = toc;
horzcat('This script took ', num2str(time),' to run.')