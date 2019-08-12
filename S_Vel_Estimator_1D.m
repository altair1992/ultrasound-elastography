function [gamma_0_1, C] = S_Vel_Estimator_1D(X,U,V)
% VELOCITYESTIMATOR Estimates velocity from an RF-signal
%   [gamma_0_1, C] = VELOCITYESTIMATOR(X, Constants, Parameters) returns
%   a matrix gamma_0_1 containing axial velocity estimates between two 
%   consecutive axial samples. C is a matrix containing the correlation
%   coefficient for the returned samples. 
%
%   PARAMETERS is a struct of type VelocityEstimatorParams.
%   CONSTANTS  is a struct of type ImagingParams. 
%
%   twoDimEn        = If set to 1, enables 2-dimensional velcity estimator
%   unwrapEn        = If set to 1, enables phase unwrapping
%   fDemodulation   = used only for 1-dimension velocity estimator,
%                     this is the fixed demodulation frequency
%   rangleGate(U)   = # of lateral samples used in each estimate
%   lateralGate(V)  = # of depth samples used in each estimate
%   ensembleLength  = # of time samples used in each estimate
%

twoDimEn = 0;
unwrapEn = 1;
fDem     = 0.5;

c = 1540;
fs =28e6;

% U = 75;
% V = 4;
O = 2;

X_conj = conj(X);

gamma_0_1 = sum(X(1:end-1,:,1:end-1).*X_conj(1:end-1,:,2:end),3);
gamma_0_1 = conv2(gamma_0_1, ones(U,1), 'valid');
gamma_0_1 = conv2(gamma_0_1, ones(1,V), 'valid');

if twoDimEn
    gamma_1_0 = sum(X(1:end-1,:,:).*X_conj(2:end,:,:),3);
    gamma_1_0 = conv2(gamma_1_0, ones(U,1), 'valid')/U;
    gamma_1_0 = conv2(gamma_1_0, ones(1,V), 'valid')/V;
end

C = sum(X(1:end-1,:,:).*X_conj(1:end-1,:,:),3);
C = conv2(C, ones(U,1), 'valid');
C = conv2(C, ones(1,V), 'valid');
C = (O/(O-1))*abs(gamma_0_1)./C;

if unwrapEn
    angle_gamma_0_1=unwrap(angle(gamma_0_1));
else
    angle_gamma_0_1=(angle(gamma_0_1));
end

% % Compute 2D Autocorrelation
% gamma_0_1= angle_gamma_0_1./(angle(gamma_1_0));
% Compute 1D Autocorrelation
gamma_0_1= angle_gamma_0_1;
% if twoDimEn
%     gamma_0_1= angle_gamma_0_1./(angle(gamma_1_0));
% else
%     gamma_0_1= 0.5*c*angle_gamma_0_1./(2*pi*fDem*fs);
% end
