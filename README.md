# Statistical Mixture Model I
## Introduction
This repository contains MATLAB codes for Inverse-Gamma Normal mixture model.  
Other codes for mixture model application for supercritical fluids (Lognormal-Nakagami mixture, Probabilistic classification) will be uploaded soon.  
If you use these codes for the publication, please cite one of the following articles.  
[1] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y.-W. (2017),The Journal of Supercritical Fluids, 119, 36-43.  
[2] Yoon, T. J., Ha, M. Y., Lee, W. B., & Lee, Y.-W. (2017), The Journal of Supercritical Fluids, 130, 364-372.  

If you have any questions, please contact me with my e-mail or researchgate account.
e-mail: oheim124@gmail.com

The specific information of these codes is as follows.

## invgampdf.m
### Descriptions
#### Input: sample values (*x*), arithmetic mean (*mu*) and standard deviation (*sigma*)
#### Output: likelihood of Inverse Gamma distribution (*y*)  
This code calculates the likelihood of the Inverse Gamma distribution. Please refer to [Wikipedia page](https://en.wikipedia.org/wiki/Inverse-gamma_distribution "Wikipedia") for the description of Inverse Gamma distribution.Population parameters of the Inverse Gamma distribution are originally shape parameter and scale parameter, but this code is written with arithmetic mean and standrad deviation as the input parameters because of the practical purposes. You can modify the code based on the Wikipedia webpage.  
### Example
This codeblock draws an Inverse Gamma distribution whose mean is 5 and standard deviation is 2.  
```MATLAB
x = 0.01:0.01:20;
y = invgampdf(x,5,2);
plot(x,y)
```

## invgamfit.m
### Descriptions
#### Input: population data (*D*)
#### Output: mean (*mu*) and standard deviation (*sigma*) of the Inverse Gamma distribution.
This code estimates the population parameters of the Inverse Gamma distribution based on the sample data.
### Example
```MATLAB
data = gamrnd(10,20,[10000,1]); % random number generation sampled from the Gamma distribution
data = data.^(-1);
m1 = mean(data);
s1 = std(data);
% So far, we generated the samples which follow the Inverse Gamma distribution whose arithmetic 
% mean and standard deviation is m1 and s1.
[mu,sigma]=invgamfit(data); % MLE for the Inverse Gamma distribution.
% You can check if the algorithm worked by comparing (m1, s1) with (mu and sigma)
% Let us check if the algorithm worked.
histogram(data,'Normalization','pdf')
hold on
x = 1E-4:1E-4:0.02;
y = invgampdf(x,mu,sigma);
plot(x,y);
```

## IGN.m
### Descriptions
#### Input: population data (*values*), covergence tolerance (*epsilon*), number of sampled data (*K*)
#### Output: arithmetic means (*mu*), standard deviations(*sigma*), and mixture weights(*p*) of the Inverse Gamma and the Normal distributions. 
This code calculates the population parameters of the Inverse Gamma-Normal mixture model using the Method of Moments.  
If the convergence rate is too slow because of a large number of population data, they can be randomly sampled.

### Example
In this example, we generate the synthetic data which consist of the random numbers sampled from the Inverse Gamma distribution and the Normal distribution.
#### Step 1: Generate random numbers.
```MATLAB
data1 = gamrnd(10,20,[100000,1]); % random numbers sampled from the Gamma distribution
data1 = data1.^(-1); % Mean and standard deviation of these random numbers are 0.0055 and 0.0020, respectively.
data2 = normrnd(0.012,0.002,[100000,1]); % random numbers sampled from the Normal distribution
data = [data1;data2];
histogram(data,'Normalization','pdf')
hold on
```
#### Step 2: Apply the Inverse Gamma-Normal mixture model.
We apply the Method of Moments to the data.
```MATLAB
[mu,sigma,p,counter]=IGN(data,1E-7,100000);
x = 1E-4:1E-4:.02; 
y1 = p(1)*invgampdf(x,mu(1),sigma(1));
y2 = p(2)*normpdf(x,mu(2),sigma(2));
plot(x,y1); plot(x,y2); plot(x,y1+y2);
```
By comparing the obtained parameters and probability distribution functions with the original data, you can know that the mixture model was well fitted to the synthetic data.
## IGN2.m
### Descriptions
#### Input: population data (*values*), covergence tolerance (*epsilon*), number of sampled data (*K*)
#### Output: arithmetic means (*mu*), standard deviations(*sigma*), and mixture weights(*p*) of the Inverse Gamma and the Normal distributions. 
This code calculates the population parameters of the Inverse Gamma-Normal mixture model using the Method of Moments.  
If the convergence rate is too slow because of a large number of population data, they can be randomly sampled.

### Example
The synthetic data obtained in the previous example is used.
#### Apply the Inverse Gamma-Normal mixture model.
We apply the Maximum Likelihood Estimation (MLE) to the data.
```MATLAB
data1 = gamrnd(10,20,[100000,1]); % random numbers sampled from the Gamma distribution
data1 = data1.^(-1); % Mean and standard deviation of these random numbers are 0.0055 and 0.0020, respectively.
data2 = normrnd(0.012,0.002,[100000,1]); % random numbers sampled from the Normal distribution
data = [data1;data2];
histogram(data,'Normalization','pdf')
hold on
[mu, sigma, p, counter, loglikelihood]=IGNMLE(data,1E-7,10000);
x=1E-4:1E-4:.02;
>> y1=p(1)*invgampdf(x,mu(1),sigma(1));
y2=p(2)*normpdf(x,mu(2),sigma(2));
plot(x,y1)
hold on
plot(x,y2)
plot(x,y1+y2)
```
By comparing the obtained parameters and probability distribution functions with the original data, you can know that the mixture model was well fitted to the synthetic data.
