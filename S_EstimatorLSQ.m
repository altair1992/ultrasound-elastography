function [S, S_C] = S_EstimatorLSQ(V, C, delta_m)
%S_EstimatorLSQ Least-squares implementation of strainEstimator (Simple version) 
%
%  
%
%   See also strainEstimator, VELOCITYESTIMATOR, StrainEstimatorParams, ImagingParams 

%delta_m = 125;               % For 20MHz 90, 28MHz 125, 60MHz 270, 100MHz 450

kappa      = 1;

n = delta_m+1;

h   = 12/(n*(n^2-1)) * ( (n:-1:1) - (n+1)/2 )';
h_c = 1/n*ones(n,1);

S   = kappa*conv2(V,h,'valid');
S_C = conv2(C,h_c,'valid');
