function theta = project3_2(theta,mu,dy,v_0,X,N,max_iter)

k = 1;
err = zeros(max_iter,1);
flag = 0;

warning('off','all');

while(k <= max_iter)
    theta_k = theta;
    y = -p13J(theta,mu,dy,v_0,X,N)\p13F(theta,mu,dy,v_0,X,N);
    theta = theta + y;
    err(k) = norm(theta_k - theta,inf);
    if(err(k) < 1e-8)
        flag = 1;
        break
    end
    k = k + 1;
end


if(flag == 1)
    fprintf('k      ||theta_k - theta_k-1||_inf\n');
    for j = 1:k
        fprintf('%d           %12.9f\n',j+1,err(j));
    end
    x = [0; cumsum(dy*tan(theta))]; 
    y = -dy * (0:N)'; 
    plot(x,y,'.-');
    axis equal;
    drawnow
else 
    fprintf('No Convergence\n')
end

end
