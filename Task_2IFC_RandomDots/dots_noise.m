% This code reproduces the analyses in the paper
% Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven 
% by decision uncertainty and alters serial choice bias. 
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai, 2016
% anne.urai@gmail.com

function [ stimuli] = dots_noise(dots, NFR)

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

% preallocate
stimuli = nan(NFR, 2, NDOTS);

% generate random starting points within a circular aperture
rad     = RADIUS*sqrt(rand(NDOTS,1));
theta   = 2*pi*rand(NDOTS,1);
pos     = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates

for frame = 1:NFR,
    
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
    stimuli(frame, :, :)  = pos';

end

end