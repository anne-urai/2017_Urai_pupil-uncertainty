function [ window] = dots_noise_draw(window, dots)

% Generates a noisy dot patch according to a limited lifetime
% algorithm, see:
% Pilly, P. K., & Seitz, A. R. (2009).
% What a difference a parameter makes: a psychophysical comparison
% of random dot motion algorithms. Vision Research, 49(13), 1599?612.
%--------------------------------------------------------------------------

NDOTS       = dots.nDots; % 774, 3097
RADIUS      = dots.radius;
INNER       = dots.innerspace; %inner circle doesn't contain dots
COH         = 0;

% generate random starting points within a circular aperture
rad     = RADIUS*sqrt(rand(NDOTS,1));
theta   = 2*pi*rand(NDOTS,1);
pos     = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates

% define the noise dots
noisedots           = round(NDOTS*(1-COH)); %nr of noisedots based on the coherence
t1                  = rand(NDOTS,1);
[t1,t2]             = sort(t1);
noiseindex          = t2(1:noisedots); %random subset of dots

% replot the noisedots somewhere in the aperture
rad                 = RADIUS*sqrt(rand(noisedots,1));
theta               = 2*pi*rand(noisedots,1);
[pos(noiseindex, 1), pos(noiseindex,2)] = pol2cart(theta, rad);

%find the dots that have left the aperture
outindex            = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) >= RADIUS); %index dots that are outside the aperture
pos(outindex, 2)    = -pos(outindex,2); %move back to the top, only change y coordinate

%find the dots that are too close to the fixation (1 degree)
innerindex          = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) <= INNER);
rad                 = RADIUS*sqrt(rand(length(innerindex),1));
theta               = 2*pi*rand(length(innerindex),1);
[pos(innerindex, 1), pos(innerindex,2)] = pol2cart(theta, rad);

% do this again to make sure there are none left
innerindex          = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) <= INNER);
rad                 = RADIUS*sqrt(rand(length(innerindex),1));
theta               = 2*pi*rand(length(innerindex),1);
[pos(innerindex, 1), pos(innerindex,2)] = pol2cart(theta, rad);

% save
stimuli(:, :)  = pos';

% draws the dots on the screen for this frame
Screen('DrawDots',window.h,squeeze(stimuli( :, :)), ...
    dots.size, dots.color, window.center, 2);
end

