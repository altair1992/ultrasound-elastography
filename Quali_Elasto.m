function [SNRe, CNRe] = Quali_Elasto(S_Frame,pos_target, pos_background,scan_mode)
x_target = pos_target(1); y_target = pos_target(2);
x_back = pos_background(1); y_back = pos_background(2);

A = S_Frame;
% Compute Signal-to-Noise for Elastogram
m = mean(A(:));
sig = std(A(:));

SNRe = m/sig;

switch (scan_mode)
    case 'curve'
        % Compute Contrast-to-Noise for Elastogram (Data Set 1)
        x_ws = 30;              % Comp. exam = 75
        y_ws = 4;               % Comp. exam = 2
    case 'linear'
        % For Data Set 2
        x_ws = 24;              
        y_ws = 3;
    case 'virtual'
        % Compute Contrast-to-Noise for Elastogram (Virtual Phantom)
        x_ws = 24;              
        y_ws = 3;
end

B = A(x_back-x_ws: x_back+x_ws, y_back-y_ws: y_back+y_ws);
T = A(x_target-x_ws: x_target+x_ws, y_target-y_ws: y_target+y_ws);

m_back = mean(B(:)); m_targ = mean(T(:));
var_back = (std(B(:)))^2; var_targ = (std(B(:)))^2;

CNRe = 2*(m_back - m_targ)^2/ (var_back + var_targ);

end