function [elastodata, alfa] = make_elastography(dtau,threshold,nframes,smoothweight,kernel,order)
% function [elastodata, alfa] = make_elastography(dtau,threshold,nframes,smoothweight,kernel,order)
%
% Computes thresholded and smoothed outputs of the strain data which gives
% an elastogram. Thresholding is done on individual frames to maximize the
% dynamic area of each frame. Smoothing is performed both in the time and
% the spatial domain to suppress noise.
%
% Inputs
% dtau: 3D matrix of strain values
% threshold: 1x2 matrix of initial threshold values
% nframes: Number of frames for time domain smoothing
% weight: Weight of frames for time domain smoothing, gauss or average
% kernel: Mask for spatial smoothing, uniform, gauss, kalman or median
% order: Size of the kernel matrix
%
% Outputs
% elastodata: 3D matrix of elastography data
% alfa: Threshold values for individual frames


if nargin < 6
    order = 5;              % Default 5
end

if nargin < 5
    kernel = 'median';
end

if nargin < 4
    smoothweight = 'average';
end

if nargin < 3
    nframes = 5;
end

if nargin < 2
    threshold = [-5e-3 5e-3];
end

elastodata = dtau;
[samples,lines,frames,iterations] = size(elastodata);

%Change the dynamic area of the data by thresholding values
min_threshold = threshold(1);
max_threshold = threshold(2);

%Initially threshold entire dataset before individual thresholding of frames
elastodata(elastodata > max_threshold) = max_threshold;
elastodata(elastodata < min_threshold) = min_threshold;


%Number of bars for computing individual histogram
bars = 400;

%Alpha-low and alpha-high for each frame is saved in variable alfa
alfa = zeros(frames,iterations,2);

%Individual thresholding of frames
for frame = 1:frames

    for iteration = 1:iterations

        frame1 = elastodata(:,:,frame,iteration);

        %Calculate histogram of individual frames
        [N x_out] = hist(frame1(:),bars);
        N(N==max(N)) = max_threshold;

        %Computing new histogram leaving out already thresholded values
        N = N(2:end-1);x_out=x_out(2:end-1);

        %Find max value and standard deviation
        maxVal = max(N);
        stdVal = std(N);

        %Define threshold value as max value plus/minus 3*standard deviation
        tVal = maxVal-3*stdVal;

        %Find all values above threshold value
        aboveThreshold = find(N>tVal);

        %Define alfaLow and alfaHigh
        alfaLow = x_out(aboveThreshold(1));
        alfaHigh = x_out(aboveThreshold(end));

        %Threshold values that are above or below alpha values
        frame1(frame1 < alfaLow) = alfaLow;
        frame1(frame1 > alfaHigh) = alfaHigh;

        %Replace original frame with thresholded frame
        elastodata(:,:,frame,iteration) = frame1;

        %Save alpha values
        alfa(frame,iteration,:) = [alfaLow alfaHigh];

    end
end

%Computing filter mask
if strcmp(kernel,'gauss')
    zfactor = round(samples/lines);
    mask = gausswin(order*zfactor)*gausswin(order)';
    mask = mask/sum(mask(:));

elseif strcmp(kernel,'uniform')
    zfactor = round(samples/lines);
    mask = ones(order*zfactor,1)*ones(1,order);
    mask = mask/sum(mask(:));

elseif strcmp(kernel,'median')
    zfactor = round(samples/lines);
    medsize = [order*zfactor order];
    
elseif strcmp(kernel,'kalman')
    zfactor = round(samples/lines);
    medsize = [order*zfactor order];
    
else
    error('Use gauss, uniform or median for kernel')
end

%Compute middle value of number of frames
middleVal = floor(nframes/2);

%Make time-smoothing statement to be executed based on numFramesSmooth
timeSmooth = 'img = ';
for l=1:nframes
    timeSmooth = [timeSmooth 'weight(' num2str(l)...
        ')*elastodata(:,:,k-(' num2str(middleVal-l+1) '))+'];
end
timeSmooth = [timeSmooth(1:end-1) ';'];

%Compute gaussian or average weight for frames
if strcmp(smoothweight,'gauss')
    weight = gausswin(nframes)/sum(gausswin(nframes));
elseif strcmp(smoothweight,'average')
    weight = ones(nframes)/sum(ones(nframes));
else
    error('Use gauss or average for smoothing weight of frames')
end

%Persistence value for frames at the very beginning or very end of sequence
perVal = 0.35;
tic
for k=1:frames

    for iteration=1:iterations

        %Assign value to first frame of sequence
        if k==1
            img = elastodata(:,:,k,iteration);

        %Time smoothing over several frams
        elseif k > middleVal
            %Frames at very end of sequence
            if k > frames-middleVal
                img_k_1 = elastodata(:,:,k-1,iteration);
                img = elastodata(:,:,k,iteration);
                %img = (1-perVal)*img + perVal*img_k_1;

            %Execute time smoothing statement for all other frames
            else
                eval(timeSmooth)
            end

            %Frames at very beginning of sequence
        else
            img_k_1 = elastodata(:,:,k-1,iteration);
            img = (1-perVal)*img + perVal*img_k_1;
        end

        %Spatial smoothing
        if strcmp(kernel,'median')
            img = medfilt2(img,medsize);
%             U = floor(medsize(1)*0.5);
%             V = 5;
%             img = conv2(img,ones(U,1),'same')/U;
%             img = conv2(img,ones(1,V),'same')/V;
        elseif strcmp(kernel,'kalman')
            %img = img;
            img = kalman(img);
            %img = medfilt2(img,medsize);
        else
            img = conv2(img,mask,'same');
        end

        %Store result
        [m,n]=size(img);
        elastodata(1:m,1:n,k,iteration) = img;

    end
end
toc