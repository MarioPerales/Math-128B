function [x, w] = nonlinear_shooting(f, dfdy, dfdyp, a, b, alpha, beta, N, tol, M)
%NONLINEAR_SHOOTING From section 11.2 of Burden/Faires, as system of 4 eqns
%   [X, W] = NONLINEAR_SHOOTING(F, DFDY, DFDYP, A, B, ALPHA, BETA, N, TOL, M)

h = (b-a) / N;
tk = (beta - alpha) / (b-a);
x = a + h*(0:N);

fprintf(' k   |w_{1,n} - beta|   \n');
fprintf('------------------------\n');
for k = 1:M
    w = [[alpha; tk; 0; 1], zeros(4,N)];
    for i = 1:N
        F = @(x, u) [u(2); f(x, u(1), u(2)); u(4); ...
                     dfdy(x, u(1), u(2))*u(3) + dfdyp(x, u(1), u(2))*u(4)];
        k1 = h*F(x(i), w(:,i));
        k2 = h*F(x(i) + h/2, w(:,i) + k1/2);
        k3 = h*F(x(i) + h/2, w(:,i) + k2/2);
        k4 = h*F(x(i) + h, w(:,i) + k3);
        w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
    fprintf('%2d     %12.8f\n', k, abs(w(1,end) - beta));
    if abs(w(1,end) - beta) < tol, return; end
    if k == M, error('No convergence.'); end
    tk = tk - (w(1,end) - beta) / w(3,end);
end
