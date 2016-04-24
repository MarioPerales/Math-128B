tic;
n = 50;
k = 1000;
omega = 2/(1+sqrt(8)/(n+1));

A = delsq(numgrid('S', n+2)) * (n+1)^2;
b = ones(n^2,1);

[~, jacobi_error] = jacobi_method(A,b,k);
[~, gauss_error] = gauss_method(A,b,k);
[~, SOR_error] = SOR_method(A,b,omega,k);
[~, cg_error] = my_cg(A,b,k);
[~, pre_cg_error] = pre_cg(A,b,k);
[~, pre_cg_ichol_error] = pre_cg_ichol(A,b,k);

semilogy(jacobi_error(1,:),'r');
hold on
semilogy(gauss_error(1,:),'g');
hold on
semilogy(SOR_error(1,:),'b');
hold on 
semilogy(cg_error(1,:), 'k');
hold on 
semilogy(pre_cg_error(1,:), 'm');
hold on
semilogy(pre_cg_ichol_error(1,:), 'c');
hold off
legend('Jacobi','Gauss-Seidel','SOR','CG', 'Pre-CG', 'Pre-CG (with ichol function)','Location','east');
ylabel('Relative Error');
xlabel('Number of iterations');

time = toc;
horzcat('This script took ', num2str(time), ' seconds to run.')