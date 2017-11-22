function [mu, sigma, p, counter, loglikelihood]=ignmle(values,epsilon,K)
% This function is used to fit the inverse gamma distribution to data.
% This code can be freely distributed.
% Please cite one of the following article if you use employ this code.
% [1] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y. W. (2017). 
%     The Journal of Supercritical Fluids, 119, 36-43.
% [2] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y. W. (2017).
%     The Journal of Supercritical Fluids, 130, 364-372.

% random sampling K values
values=datasample(values,K);
% values should be nX1 vector.
% initialize
idx = kmeans(values,2);
if mean(values(idx==1)) > mean(values(idx==2))
    mu(1) = mean(values(idx==2));
    mu(2) = mean(values(idx==1));
    sigma(1) = std(values(idx==2));
    sigma(2) = std(values(idx==1));
else
    mu(1) = mean(values(idx==1));
    mu(2) = mean(values(idx==2));
    sigma(1) = std(values(idx==1));
    sigma(2) = std(values(idx==2));
end
values=values';
counter = 0;
difference = epsilon;
p = ones(2,1)/2;

% Initial guess
alpha = (mu(1)/sigma(1))^2+2;
if alpha > 60
    alpha = 60;
    mu(1) = sqrt(alpha-2)*sigma(1);
end
    
% now iterate.
while (difference >= epsilon)
    % E step: soft classification of the values into one of the mixtures
    class(1,:) = p(1)*invgampdf(values,mu(1),sigma(1));
    class(2,:) = p(2)*normpdf(values, mu(2), sigma(2));
    
    % Normalize.
    class = class./repmat(sum(class),2,1);
    
    % M step: Maximum Likelihood estimate the parameters of each class
    % (i.e., p, mu, sigma)
    mu_old = mu;
    sigma_old = sigma;
    p_old=p;
    
    % inverse gamma distribution maximization
    % method of moments
    mu(1) = sum(class(1,:).*values) / sum(class(1,:));
    sigma(1) = sqrt( sum(class(1,:).*(values - mu_old(1)).^2) /  sum(class(1,:)) );
    alpha_old = (mu(1)/sigma(1))^2+2;
                
    % iteration to obtain beta and alpha   
    for j=1:100
        C = -log(sum(class(1,:)./values))-sum(class(1,:).*log(values))/sum(class(1,:));
        alpha = invpsi(log(sum(class(1,:))*alpha_old)+C);
        if alpha_old-alpha < 10^(-3)
            break;
        end
        alpha_old=alpha;
    end    
    beta = alpha*sum(class(1,:))/sum(class(1,:)./values);
    mu(1) = beta/(alpha-1);
    sigma(1) = beta/((alpha-1)*sqrt(alpha-2));
    p(1) = mean(class(1,:));
    
    % normal distribution maximization
    mu(2) = sum(class(2,:).*values) / sum(class(2,:));
    sigma(2) = sqrt( sum(class(2,:).*(values - mu(2)).^2) /  sum(class(2,:)) );
    p(2) = mean(class(2,:));
    
    difference(counter+1) = sum(abs(mu_old - mu)) + sum(abs(sigma_old - sigma)) + ...
      sum(abs(p_old - p));
    counter = counter +1;    
    loglikelihood = sum(class(1,:).*log(p(1)*invgampdf(values,mu(1),sigma(1))))+sum(class(2,:).*log(p(2)*normpdf(values,mu(2),sigma(2))));
    if p(1) < 2E-2 || p(2) < 2E-2
        break;
    end
end

    function Y=invpsi(X)
        L = 1;
        Y = exp(X);        
        while L > 10e-6
            Y = Y + L*sign(X-psi(Y));
            L = L / 2;
        end
    end
end
