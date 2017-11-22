function [mu,sigma]=invgamfit(D)
% This function is used to fit the inverse gamma distribution to data.
% This code can be freely distributed.
% Please cite one of the following article if you use employ this code.
% [1] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y. W. (2017). 
%     The Journal of Supercritical Fluids, 119, 36-43.
% [2] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y. W. (2017).
%     The Journal of Supercritical Fluids, 130, 364-372.

mu=mean(D); sigma2=var(D);
alpha_old=mu^2/sigma2+2;
C = -log(sum(1./D))-mean(log(D));

for j=1:100        
    alpha = invpsi(log(length(D)*(alpha_old))+C);
    if alpha_old-alpha < 10^(-5)
        break;
    end
    alpha_old=alpha;
end

beta = alpha*length(D)/sum(1./D);
mu = beta/(alpha-1);
sigma2 = (beta/((alpha-1)*sqrt(alpha-2)))^2;
sigma = sqrt(sigma2);

    function Y=invpsi(X)
        L = 1;
        Y = exp(X);        
        while L > 10e-6
            Y = Y + L*sign(X-psi(Y));
            L = L / 2;
        end
    end
end
