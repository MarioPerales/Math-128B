function theta = project3_1(theta,mu,dy,v_0,X,N,max_iter)

k = 1;

while(k <= max_iter)
    theta_k = theta;
    y = -jacobian_13(theta,mu,dy,v_0,N)\F_13(theta,mu,dy,X,v_0,N);
    theta = theta + y;
    k = k + 1;
    if(norm(theta_k - theta,inf) <= 1e-2)
        break
    end   
end

