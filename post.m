% Post_process
function [SNRe_i, CNRe_i]= post(pos_target,pos_back,scan_mode,fs1)
%path(path, './Results');

load('./Results/elastodata.mat')
load('./Results/velocity.mat')
load('./Results/eyy.mat')

% Set up the result movie
writerObj = VideoWriter('./Results/out.avi');
writerObj.FrameRate = 2;
open(writerObj);

% 
iterations = 1;
%scan_mode = 'curve';             % Scan mode = 'linear', 'curve'
%scan_mode = 'linear';             % Scan mode = 'linear', 'curve'
kappa = (1540/2)*100*(1/20000000);
trans_pitch = 3e-4 *100;

%fs1 = 28e6;
kappa1 = (1540/2)*100*(1/fs1);
trans_pitch1 = 5.133e-4 *100;

for idx = 1:size(elastodata,3)
    elasto(:,:,idx) = elastodata(:,:,idx,iterations);
    figure(1)
    frame = elasto(:,:,idx);
    [SNRe, CNRe] = Quali_Elasto(frame,pos_target, pos_back,scan_mode);
    SNRe_i(idx) = SNRe; CNRe_i(idx) = CNRe;
    switch (scan_mode)
        case 'linear'
            %imagesc(frame)
            X = trans_pitch*(1:size(elastodata,2));
            Y = kappa*(1:size(elastodata,1));
            pcolor(X,Y,frame)
            set(gca,'YDIR','reverse')           % Reverse Y-Axis
            rectangle('Position',[1.26 1.06 0.18 0.18])   % For Target area
            rectangle('Position',[1.26 0.1025 0.18 0.18]) % For Background area
        case 'curve'
            [X, Y, data] = UWscanconv_result(frame,1,-15.5,31,4,1.4);
            pcolor(X,Y,data)
%             rectangle('Position',[-0.2 -1.1 0.4 0.4])
%             rectangle('Position',[-0.2 -2.8 0.4 0.4])
        case 'virtual'
            X = trans_pitch1*(1:size(elastodata,2));
            Y = kappa1*(1:size(elastodata,1));
            pcolor(X,Y,frame)
            set(gca,'YDIR','reverse')           % Reverse Y-Axis
            rectangle('Position',[2.0 1.8 0.18 0.18])   % For Target area
            rectangle('Position',[2.0 0.8 0.18 0.18]) % For Background area
            
    end
    xlabel('Width [cm]')
    ylabel('Depth [cm]')
    %caxis([-35,35])
    title(['Elastogram (idx:', num2str(idx),',SNRe:',num2str(SNRe),',CNRe:',num2str(CNRe),')'])
    shading('interp')    
    colormap(jet)
    axis('equal')
    colorbar('location','EastOutside')
    
    frame = getframe(gcf);
    writeVideo(writerObj,frame);    
    
    % Save image
    if idx < 10
        fignum = strcat('000',num2str(idx));
    elseif idx < 100
        fignum = strcat('00',num2str(idx));
    else idx < 1000
        fignum = strcat('0',num2str(idx));
    end
    %print('-f1','-djpeg','-r600', sprintf('./Results/Fig%4s',fignum));    % Save the figure to file
    close figure 1  
end
close(writerObj)
% Save Velocity and Strain Frame
idx_show = 3;
vel_show = Vel_frame(:,:,idx_show);
eyy_show = S_frame(:,:,idx_show);

figure(2)
switch (scan_mode)
    case 'linear'
        X = trans_pitch*(1:size(vel_show,2));
        Y = kappa*(1:size(vel_show,1));
        pcolor(X,Y,vel_show)
    case 'curve'
        [X, Y, data] = UWscanconv_result(vel_show,1,-15.5,31,4,1.4);
        pcolor(X,-Y,data)
    case 'virtual'
        X = trans_pitch1*(1:size(vel_show,2));
        Y = kappa1*(1:size(vel_show,1));
        pcolor(X,Y,vel_show)
            
end

% [cMin, cMax] = autoScale(disp_show, 'vel');
% caxis([cMin, cMax]); 

title(['Velocity (idx:', num2str(idx_show),')'])
ylabel('Depth [cm]')   
xlabel('Width [cm]') 
colormap(jet)
axis('equal')
shading('interp')
colorbar('location','EastOutside')
set(gca,'YDIR','reverse')           % Reverse Y-Axis
%print('-f2','-djpeg','-r600', sprintf('./Results/Fig_vel%4s',num2str(idx_show)));

figure(3)
switch (scan_mode)
    case 'linear'
        X = trans_pitch*(1:size(eyy_show,2));
        Y = kappa*(1:size(eyy_show,1));
        pcolor(X,Y,eyy_show)
    case 'curve'
        [X, Y, data] = UWscanconv_result(eyy_show,1,-15.5,31,4,1.4);
        pcolor(X,-Y,data)
    case 'virtual'
        X = trans_pitch1*(1:size(eyy_show,2));
        Y = kappa1*(1:size(eyy_show,1));
        pcolor(X,Y,eyy_show)
end

[cMin, cMax] = autoScale(eyy_show, 'strain');
caxis([cMin, cMax]); 
%

%caxis([0, cMax]); 
title(['Axial Strain (idx:', num2str(idx_show),')'])

ylabel('Depth [cm]')   
xlabel('Width [cm]') 
colormap(jet)
axis('equal')
shading('interp')
colorbar('location','EastOutside')
set(gca,'YDIR','reverse')           % Reverse Y-Axis
%print('-f3','-djpeg','-r600', sprintf('./Results/Fig_eyy%4s',num2str(idx_show)));

end
