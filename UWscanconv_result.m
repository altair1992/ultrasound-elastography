function [X, Y, data] = UWscanconv_result(data, ArrayType, AcquiredLateralMin, AcquiredLateralSpan, ApexVerticalCm, ApexLateralCm, DecimationFactor)
%ANTARESSCANCONVERT(DATA, ArrayType, AcquiredLateralMin, AcquiredLateralSpan, ApexVerticalCm, ApexLateralCm)
% Displays the scan-converted image of DATA, depending on the type of array transducer used.
% Currently handles linear, curved, and phased arrays.
% Currently handles parallelogram (including rectangular) and sector views.  
% This program does not generate the envelope of RF, does not apply colorflow or power colormap.  
% It uses DATA as it is (intesity-wise) for display. 
%
% Inputs:
% DATA is whatever 2D image you want to scan-convert
% ArrayType:   0 Linear Array (parallelogram)
%              1 Curved Array (sector)
%              2 Flat Phased Array (sector)
%              3 Linear-phased array (sector)
% AcquiredLateralMin: 
%      for curved and phased array: angle of first vector w/ vertical axis, in degrees.
%      for linear array: steering angle, in radian
% AcquiredLateralSpan: from DUI (in cm for linear arrary, degrees for phased arrays)
% ApexVerticalCm: from DUI
% ApexLateralCm: from DUI
%	manually recorded from DUI:
%		sys -> param -> hardwarecontroll -> hardwarealgorith -> cbimages -> cbimages0 -> scalaroptions
% (AcquiredAxialSpanCm not used because depth is calculated based on 1540m/s and 40MHz sampling rate.)
% DecimationFactor: decimation factor used on I/Q for flow processing
%
%Note:  This program assumes the data collect starts at depth=0 (element surface).
% So in DUI: main -> sys -> be -> rfdata options -> start 0.000000
%
% This is part of the URI offline processing tools, ver 0.6x
% See also URIBSCANCONV, URIFULLFOVCOLOR, URIFULLFOVCOLOR, URICOLORMAP, URIPOWERMAP

% Jerome Mai, jjmai@ucdavis.edu, 12/11/02
% copyright UC Davis

f_mode = 0;     % figure mode (On = 1, Off = 0)

if ~exist('DecimationFactor', 'var')
   DecimationFactor = 1;
elseif DecimationFactor==0
   DecimationFactor = 1;
end


Nsamples = size(data,1); %number of samples in each vector
Nvectors = size(data,2); %number of vectors in the data set
depth = 5.2;             %depth of data from Tx surface in cm
%AcquiredLateralMin = -AcquiredLateralSpan/2;	%assuming FOV is centered, pointing downward

if (AcquiredLateralMin==0)&(AcquiredLateralSpan==0)
   ArrayType = 0;
end


switch ArrayType
   
case 0	%Linear Array
   SteeringAngle = AcquiredLateralMin;	%in radian
   if ~(SteeringAngle+AcquiredLateralSpan)
      %if SteeringAngle is 0, then use imagesc, else use pcolor
      figure,  imagesc(data)	%since we don't know the size of the ROI, we won't change boxaspectratio
      %   set(gca, 'PlotBoxAspectRatio', [AcquiredLateralSpan AcquiredAxialSpan 1]);
		%   set(gca, 'PlotBoxAspectRatio', [AcquiredLateralSpan depth 1]);
      %   disp('If this image looks weird, that''s because you put input degree instead of cm');
   else
      X = (((1:Nvectors)'-1)*AcquiredLateralSpan/Nvectors*ones(1,Nsamples))' +...
         sin(SteeringAngle)*((1:Nsamples)'-1)/Nsamples*depth*ones(1,Nvectors);
      Y = cos(SteeringAngle)*((1:Nsamples)'-1)/Nsamples*depth*ones(1,Nvectors);
      if (f_mode == 1)
      figure,
          pcolor(X,Y, data)
            shading('interp')
            colormap('jet')
    % 		set(gca, 'color', 'k', 'xcolor', [.3 .3 .3], 'ycolor', [.3 .3 .3]);
    % 		shading flat
        axis equal
        axis ij
    % 		backgroundcolor = [.3 .3 .3];
    %       whitebg(gcf, color)
      end
   end
   
case 1	%Curved Array
   radius = sqrt(ApexVerticalCm^2+ ApexLateralCm^2);	%distance from apex to Tx surface in cm, sqrt(apexverticalcm^2+ApexLateralCmcm^2)
	%%%%%%%%%%%%%%%%%%readying the display grid
	R = (((1:Nsamples)'-1)/Nsamples*depth+radius)*ones(1,Nvectors);	%distance of each pixel from apex in cm
	angle = (((1:Nvectors)'-1)*AcquiredLateralSpan/Nvectors*ones(1,Nsamples))' + AcquiredLateralMin;	%angle of each pixel from vertical -pi/2 axis
	Y = R.*sin(angle/180*pi-pi/2)+abs(ApexVerticalCm);
	X = R.*cos(angle/180*pi-pi/2);
    if (f_mode == 1)
       figure,
        pcolor(X,-Y, data)
        shading('interp')
        colormap('jet')
    %    set(gca, 'color', 'k', 'xcolor', [.3 .3 .3], 'ycolor', [.3 .3 .3]);
    % 	shading flat
       axis equal
       axis ij
        %backgroundcolor = [.3 .3 .3];
       %whitebg(gcf, color)
    end
   
case {2,3}	%phased array or linear-phased array
   %assume apex location is always relative to same probe origin, and probe origin is at 90 degree full FOV.
	angle = (((1:Nvectors)'-1)*AcquiredLateralSpan/Nvectors*ones(1,Nsamples))' + AcquiredLateralMin;	%angle of each pixel from vertical -pi/2 axis
   R = (((1:Nsamples)'-1)/Nsamples*depth)*ones(1,Nvectors) + abs(ApexVerticalCm)./cos(angle*pi/180);
	Y = R.*sin(angle/180*pi-pi/2)+abs(ApexVerticalCm);
	X = R.*cos(angle/180*pi-pi/2);
    if (f_mode == 1)
       figure,
        pcolor(X,-Y, data)
        set(gca, 'color', 'k', 'xcolor', [.3 .3 .3], 'ycolor', [.3 .3 .3]);
        shading flat
       axis equal
       axis ij
        %backgroundcolor = [.3 .3 .3];
       %whitebg(gcf, color)
    end
   
otherwise
   figure, 
   imagesc(data)
   disp('Unknown array type')

end

% colormap(gray)
colorbar
