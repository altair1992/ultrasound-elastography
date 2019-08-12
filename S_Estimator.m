function [S, S_C] = S_EstimatorLSQ(V, C, delta_m)
%  Centeral Differenc Strain Estimation for RPSE
%

%delta_m = 125;               % For 20MHz 90, 28MHz 125, 60MHz 270, 100MHz 450

kappa      = 1;

n = delta_m+1;

h   = zeros(delta_m + 1,1);
h(1) = 1/delta_m;
h(delta_m + 1) = -1/delta_m;

h_c   = zeros(delta_m + 1,1);
h_c(1) = 1/2;
h_c(delta_m + 1) = 1/2;

S   = kappa*conv2(V,h,'valid');
S_C = conv2(C,h_c,'valid');
