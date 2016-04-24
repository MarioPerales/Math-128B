K=assemble(20);
       [V,D]=eigs(K,24,'sa',struct('disp',0));
%        for ii=1:size(V,2)
%          disp(sprintf('lambda_%d = %g',ii,D(ii,ii)))
%          qdanim(V(:,ii))
%        end
       
d = diag(D);

lambda_1 = p3_code(K,d(18) - (d(19)-d(18))*.6*d(18));
lambda_2 = p3_code(K,d(18) + (d(19)-d(18))*.9*d(18));
c = (d(18) + d(19))/2;
lambda_3 = p3_code(K,c + .15*(abs(d(18) - d(19)))/2);
lambda_4 = p3_code(K,d(22) - .005*d(22));

error_1 = abs(lambda_1 - d(18));
error_2 = abs(lambda_2 - d(18));
error_3 = abs(lambda_3 - d(19));
error_4 = abs(lambda_4 - d(22));

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
title('Rayleigh Quotient Iteration');
