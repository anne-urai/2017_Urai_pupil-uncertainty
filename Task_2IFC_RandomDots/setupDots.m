function [dots, fix] = setupDots(window, setup)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOTS CHARACTERISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dots                   = struct; %preallocate
% parameters that stay constant
dots.radius            = deg2pix(window, 12); %keep constant
dots.innerspace        = deg2pix(window, 2);
dots.lifetime          = 4; % in frames
dots.nvar              = 3; %interleave two variants of the stimulus
dots.color             = [255 255 255]; % always 100% dot contrast

% determine the appearance of the dots - as in Siegel, 2007
dots.speed             = deg2pix(window, 11.5); % speed of the dots in degrees/second
dots.size              = deg2pix(window, 0.2); %size of each dot in degrees
dots.density           = 1.7; % dot density in dots per degree^2
dots.nDots             = round(dots.density*pi*pix2deg(window, dots.radius)^2); %number of dots in the circle, calculated from density (can also just be a fixed nr, eg. 500

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is the fixation descibed in Thaler et al. 2013
fix.dotsize         = deg2pix(window,0.2); % fixation dot
fix.circlesize      = deg2pix(window,0.6); % circle around fixation dot
fix.color           = [255 0 0];

end