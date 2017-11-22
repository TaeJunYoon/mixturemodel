function y=invgampdf(x,mu,sigma)
% This code can be freely distributed.
% Please cite one of the following article if you use employ this code.
% [1] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y. W. (2017). 
%     The Journal of Supercritical Fluids, 119, 36-43.
% [2] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y. W. (2017).
%     The Journal of Supercritical Fluids, 130, 364-372.
alpha = (mu/sigma)^2+2;
beta = mu*(alpha-1);
y=beta^alpha/gamma(alpha)*x.^(-alpha-1).*exp(-beta./x);
end
