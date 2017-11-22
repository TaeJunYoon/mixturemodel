function [mu,sigma,p,counter]=ign(values,epsilon,K)
% This function is used to fit the inverse gamma distribution to data.
% This code can be freely distributed.
% Please cite one of the following article if you use employ this code.
% [1] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y. W. (2017). 
%     The Journal of Supercritical Fluids, 119, 36-43.
% [2] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y. W. (2017).
%     The Journal of Supercritical Fluids, 130, 364-372.

% random sampling K values
values=datasample(values,K);
values=sort(values);
% values should be nX1 vector.
% initialize
choice = datasample(values,2);
choice = sort(choice);
mu = zeros(1,2); sigma = zeros(1,2);
for i=1:2
    mu(i)=choice(i);
    sigma(i) = std(values);
end
values=values';
counter = 0;
difference = epsilon;
p = ones(1,2)/2;
class(1,:) = p(1)*invgampdf(values,mu(1),sigma(1));
class(2,:) = p(2)*normpdf(values, mu(2), sigma(2));
class = class./repmat(sum(class),2,1);
% now iterate.
while (difference >= epsilon)    
    % M step: Maximum Likelihood estimate the parameters of each class
    % (i.e., p, mu, sigma)
    mu_old = mu;
    sigma_old = sigma;
    p_old=p;
        
    % method of moments
    for i=1:2
        mu(i)=sum(class(i,:).*values)/sum(class(i,:));
        sigma(i)=sqrt(sum(class(i,:).*(values-mu(i)).^2)/sum(class(i,:)));
    end
    
    % E step: soft classification of the values into one of the mixtures
    class(1,:) = p(1)*invgampdf(values,mu(1),sigma(1));
    class(2,:) = p(2)*normpdf(values, mu(2), sigma(2));        
    
    % Normalize.
    class = class./repmat(sum(class),2,1);
    
    for i=1:2
        p(i) = mean(class(i,:));
    end    
    difference = sum(abs(mu_old./mu-1)) + sum(abs(sigma_old./sigma-1)) + sum(abs(p_old./p-1));    
    counter = counter +1;            
end
