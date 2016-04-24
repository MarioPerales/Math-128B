function lambda = p1_code(K,mu)

[m,n] = size(K); %#ok<ASGLU>
N = size(K,1);
v0 = randn(n,1);
k = 2;
v = v0;

shifted_matrix = K - mu*speye(N,N);
w = shifted_matrix\v;
v = w/norm(w);
lambda(1) = v'*K*v;

while(k <= 10000) 
    shifted_matrix = K - mu*speye(N,N);
    w = shifted_matrix\v;
    v = w/norm(w);
%     qdplot(v);
%     drawnow;
    lambda(k) = v'*K*v; %#ok<*AGROW>
    
    if(abs(lambda(k) - lambda(k - 1)) < 1e-14)
        lambda = lambda';
        return
    end
    k = k + 1;
end