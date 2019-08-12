% Generate strain image using simple phase-based strain estimation
%                     
% Step 1: Load RF frames 
% Step 2: Compute velocity and strain field using least-square method
% Step 3: Generate results (Velocity, Strain, Elastogram)
%
% Inputs: 
%
% image_save  - 1 Save, 0 No save
%
% For data loading
%
% scan_mode(Array type) 
%             'linear'  - Load NTNU data set
%             'curve'   - Load silicon phantom test data set collecting 
%                       by a portable ultrasound device
%             'virtual' - Load a numerical phantom data set
%
% For velocity estimation
%
% U           - average kernel size for depth direction for smoothing 
%             velocity (default = 20, optimal in terms of CNR = 75 samples)
% V           - average kernel size for laternal direction for smoothing 
%             velocity (default = 4 A-lines)
%
% For Least-Square Kernel
%
% delta_m     - segment size for least-square
% linear_LS   - Least-square mode (1: 1st order, otherwise: 2nd order)
% For Post-processing
% a_scale     - Auto scale mode  (1: auto, otherwise: manual) that specifies initial threshold value
% nframes     - Number of frames for time domain smoothing
% smoothweight - Weight of frames for time domain smoothing, gauss or average
% s_kernel    - Mask for spatial smoothing, uniform, gauss, kalman or median
% s_order     - Size of the kernel matrix (default:5)
% 
% by Bruce Shin, bh3shin@uwaterloo.ca. April.2016
% copyright University of Waterloo, Canada
%
clc,clear,close all
%
% Set up parameters------------------------------------------------------

% Select data set
image_save = 0;                % image save for velocity and strain (1: save, otherwise: no)
scan_mode = 'linear';          % Array type 'linear', 'curve' 'virtual'

inc_f = 1;                      % frame increment (1 -> 0.1%, 2 -> 0.2%)
% Average kernel size for Velocity estimation
U = 20;                         % Depth direction ( U samples)
V = 4;                          % Laternal direction ( # of A-line)

% Least-Square for strain estimation  
delta_m = 125;
linear_LS = 1;                  % Least-square mode (0: Centeral Difference, 1: 1st order, 2: 2nd)

% Post-processing
a_scale = 1;                    % Auto scale mode  (1: auto, otherwise: manual) 
%                                 Specify threshold value
nframes = 2;                    % Number of frames for time domain smoothing (2 for journal paper)
smoothweight = 'gauss';         % Weight of frames for time domain smoothing, gauss or average
s_kernel = 'median';            % Mask for spatial smoothing, uniform, gauss, kalman or median
s_order = 5;                    % Size of the kernel matrix (default:5)
%--------------------------------------------------------------------------
path(path, './Data');
path(path, './Results');

% Least-Square Kernel
n = delta_m + 1;
ls_a1 = [1:n]';
if linear_LS == 1
    ls = pinv([ls_a1 ones(n,1)]);
    ls_a = -ls(1,:);        % least-square kernel for computing slope (1st order)
else
    ls =pinv([ls_a1.^2 ls_a1 ones(n,1)]);
    ls_a = -ls(2,:);       % least-square kernel for computing slope (2nd order)
end

% Step 1: Load RF Frames
%
            
switch (scan_mode)
    case 'linear'
        load('rec_data.mat')
        fs = 20e6;                          % Sampling frequency
        h = waitbar(0,'Please wait.');      % Wait bar
        frames = size(rf_data_set,3);
        % Coordinates for CNRe for Data Set 2
        pos_target = [300, 45];
        %pos_back = [900, 15];              % west side
        pos_back = [50, 45];                % north side
    case 'curve'
        load('rec_data_sonon_g.mat')          % Gelatin phantom
        %load('rec_data_sonon.mat')         % Silicon phantom
        %load('rec_data_sonon_comp.mat')         % Silicon phantom (Comprehensive exam)
        %load('rec_data_cs.mat')
        fs = 28e6;                          % Sampling frequency
        h = waitbar(0,'Please wait.');      % Wait bar
        frames = size(rf_data_set,3);
        % Coordinates for CNRe for Data Set 1
        pos_target = [1300, 33];            % For Gelatin phantom
        %pos_target = [900, 33];            % For Silicon phantom
        pos_back = [300, 33];               % north side   
        
    case 'virtual'
        load('rec_data_v_phantom_64act.mat')     % 64 active elements with FEA model at 2015.Dec
        fs = 28e6;                          % Sampling frequency
        h = waitbar(0,'Please wait.');      % Wait bar
        frames = size(rf_data_set,3);
        % Coordinates for CNRe for Data Set 2
        if fs == 100e6
            pos_target = [2500, 40];        %pos_back = [900, 15]; % west side
            pos_back = [1500, 40];          % north side
        elseif fs == 60e6
            pos_target = [1700, 40];        %pos_back = [900, 15]; % west side
            pos_back = [900, 40];           % north side
        elseif fs == 28e6
            pos_target = [750, 40];        %pos_back = [900, 15]; % west side 
            pos_back = [300, 40];           % north side
%             % When removing edge
%             pos_target = [750, 36];        %pos_back = [900, 15]; % west side 
%             pos_back = [300, 36];           % north side
        elseif fs == 40e6
            pos_target = [1000, 40];        %pos_back = [900, 15]; % west side
            pos_back = [400, 40];           % north side
        elseif fs == 50e6
            pos_target = [1400, 40];        %pos_back = [900, 15]; % west side
            pos_back = [500, 40];           % north side
        end
end    

% Step 2: Compute velocity and strain
hsize = 10; sigma = 2;
PSF = fspecial('gaussian', hsize, sigma);
tic
for frame = 1:frames-inc_f
    waitbar(frame/(frames-inc_f));
    X(:,:,1) = hilbert(rf_data_set(:,:,frame));
    X(:,:,2) = hilbert(rf_data_set(:,:,frame+inc_f));
    [Vel_frame(:,:,frame) C_Velframe(:,:,frame)] = S_Vel_Estimator(X,U,V);
    
    if linear_LS == 1
        [S_frame(:,:,frame) SC_frame(:,:,frame)]= S_EstimatorLSQ(Vel_frame(:,:,frame),C_Velframe(:,:,frame),delta_m);
    elseif linear_LS == 2
        [S_frame(:,:,frame) SC_frame(:,:,frame)]= S_EstimatorLSQ_2nd(Vel_frame(:,:,frame),C_Velframe(:,:,frame),ls_a);
    else
        % Centeral Differenc Strain Estimation for SPSE
        [S_frame(:,:,frame) SC_frame(:,:,frame)]= S_Estimator(Vel_frame(:,:,frame),C_Velframe(:,:,frame),delta_m);
    end
  
    % Build accumulated strain value
    if frame == 1
        S_A_frame(:,:,frame) = S_frame(:,:,frame);
    else
        S_A_frame(:,:,frame) = S_frame(:,:,frame-1)+S_frame(:,:,frame);
    end
end
toc
close(h)

%
% Step 3: Generate results (Velocity, Strain, Elastogram)
% 
frames = size(S_frame,3);
if image_save ==1 
    for idx=1:frames-1
        % Show the comparison result
        figure(1)
        %title('Axial displacement %s', idx)
        subplot(1,2,1)
        imagesc(Vel_frame(:,:,idx))
        xlabel('scanline')
        %caxis([-35,35])
        title(['Velocity (idx:', num2str(idx),')'])
        colormap(jet)
        colorbar('location','EastOutside')
        
        subplot(1,2,2)
        frame = S_frame(:,:,idx);
        [cMin, cMax] = autoScale(frame,'strain');
        imagesc(frame)
       % [SNRe, CNRe] = Quali_Elasto(frame(:,5:72),pos_target, pos_back,scan_mode); % remove edge
        [SNRe, CNRe] = Quali_Elasto(frame,pos_target, pos_back,scan_mode); % remove edge
        xlabel('scanline')
        title(['Strain (idx:', num2str(idx),')'])
        %title(['Strain (idx:', num2str(idx),',SNRe:',num2str(SNRe),',CNRe:',num2str(CNRe),')'])
        %caxis([-0.005,0.005])          % for instantenous strain
        %caxis([-0.02,0.02])             % for accumlated strain
        caxis([cMin, cMax]);
        colormap(jet)
        colorbar('location','EastOutside')
         
        % Save image
        if idx < 10
            fignum = strcat('000',num2str(idx));
        elseif idx < 100
            fignum = strcat('00',num2str(idx));
        else idx < 1000
            fignum = strcat('0',num2str(idx));
        end
        print('-f1','-djpeg','-r600', sprintf('./Results/Fig_v_s%4s',fignum));    % Save the figure to file
        close figure 1                
    end
end
if a_scale == 1
    [cMin, cMax] = autoScale(frame,'strain');
else
    [cMin, cMax] = [-5e-2 5e-2];
end
%
% Generate elastogram that is performed by time and spatial smoothing 
%
[elastodata alpha] = make_elastography(S_frame,[cMin cMax],nframes,smoothweight,s_kernel,s_order);
% Save results
 save('./Results/elastodata.mat','elastodata');
 save('./Results/eyy.mat','S_frame');
 save('./Results/velocity.mat','Vel_frame');
% 
 [SNRe, CNRe] = post(pos_target,pos_back,scan_mode,fs);
 save('./Results/snr_cnr.mat','SNRe','CNRe');
 
%  mean(SNRe(1:end-1))
%  mean(CNRe(1:end-1))
 