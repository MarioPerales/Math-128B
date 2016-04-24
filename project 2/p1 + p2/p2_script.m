K=assemble(20);
       [V,D]=eigs(K,24,'sa',struct('disp',0));
%        for ii=1:size(V,2)
%          disp(sprintf('lambda_%d = %g',ii,D(ii,ii)))
%          qdanim(V(:,ii))
%        end
       
d = diag(D);

mu_1 = d(18) - (d(19)-d(18))*.6*d(18); % We go slightly BELOW the eigenvalue d(18) (Why do we go LESS below than we do above? Because this eigenvalue is CLOSER to the eigenvalue below it than it is above it, so we have to compensate for this.
mu_2 = d(18) + (d(19)-d(18))*.9*d(18); % We go slightly ABOVE the eigenvalue d(18)

lambda_1 = p1_code(K,mu_1);
lambda_2 = p1_code(K,mu_2);

c = abs((d(18) + d(19)))/2;
% while(c <= d(18) + .1*d(18))
%     c = c/2;
% end

mu_3 = c + .1*(abs(d(18) - d(19)))/2; %Get it 10% from the halfway point TOWARDS the eigenvalue d(19).
mu_4 = d(22) - .005*d(22); %d(21) and d(22) are eigenvalues that are near. We go SLIGHTLY below one of them and take that as our mu.
lambda_3 = p1_code(K,mu_3);
lambda_4 = p1_code(K,mu_4);

% We compute the error based on the NEAREST eigenvalue from our mu %
error_1 = abs(lambda_1 - d(18)); 
error_2 = abs(lambda_2 - d(18));
error_3 = abs(lambda_3 - d(19));
error_4 = abs(lambda_4 - d(22));

% Plot %
semilogy(1:length(error_1),error_1,'r');
hold on
semilogy(1:length(error_2),error_2,'b');
hold on
semilogy(1:length(error_3),error_3,'g');
hold on
semilogy(1:length(error_4),error_4,'k');
legend('Part A', 'Part B', 'Part C', 'Part D');
xlabel('Number of iterations');
ylabel('Error');
title('Inverse Iteration');

% table = [mu_1, d(18); mu_2, d(18); mu_3, d(19); mu_4, d(22)];

fprintf('Mu                 |Nearest Eigenvalue        \n')
fprintf('-------------------|-------------------\n')
fprintf('%.16f |%.16f \n',mu_1,d(18))
fprintf('%.16f |%.16f \n',mu_2,d(18))
fprintf('%.16f |%.16f \n',mu_3,d(19))
fprintf('%.16f |%.16f \n\n\n\n',mu_4,d(22))

k_1 = length(error_1);
k_2 = length(error_2);
k_3 = length(error_3);
k_4 = length(error_4);

% Our actual convergence rates %
part_1_C = exp(log(error_1(k_1-1)/error_1(k_1 - 2)));
part_2_C = exp(log(error_2(k_2-1)/error_2(k_2 - 2)));
part_3_C = exp(log(error_3(k_3-1)/error_3(k_3 - 2)));
part_4_C = exp(log(error_4(k_4-1)/error_4(k_4 - 2)));

% Our estimated convergence rates %
part_1_est = (abs(mu_1 - d(18))/abs(mu_1 - d(17)))^2;
part_2_est = (abs(mu_2 - d(18))/abs(mu_2 - d(17)))^2;
part_3_est = (abs(mu_3 - d(19))/abs(mu_3 - d(18)))^2;
part_4_est = (abs(mu_4 - d(22))/abs(mu_4 - d(22)))^2;

our_data = cell(5,3);
our_data{1,1} = 'Actual Convergence';
our_data{1,2} = 'Estimated Convergence';
our_data{1,3} = 'Absolute Error';
our_data{2,1} = part_1_C;
our_data{2,2} = part_1_est;
our_data{2,3} = abs(part_1_C - part_1_est);
our_data{3,1} = part_2_C;
our_data{3,2} = part_2_est;
our_data{3,3} = abs(part_2_C - part_2_est);
our_data{4,1} = part_3_C;
our_data{4,2} = part_3_est;
our_data{4,3} = abs(part_3_C - part_3_est);
our_data{5,1} = part_4_C;
our_data{5,2} = part_4_est;
our_data{5,3} = abs(part_4_C - part_4_est);
disp(our_data)


