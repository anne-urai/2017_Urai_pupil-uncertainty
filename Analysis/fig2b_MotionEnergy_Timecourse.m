function fig2b_MotionEnergy_Timecourse()
% get the full timecourse (fix, s1, delay, s2, resp) of one trial to show
% fluctuations

if ~exist(sprintf('%s/Data/MotionEnergy/motionTimecourse.mat', mypath), 'file'),
    
    load(sprintf('%s/Data/P17/Behav/Dots_P17_s1_b1_2015-02-02_16-55-33.mat', mypath));
    
    % general experimental design
    setup.baselinecoh       = .7; %70% baseline coherence, great for MEG
    setup.trialrep          = 4; %600; %nr of trials over the whole session - 600 for MEG, more for training?
    setup.nblocks           = 1; % how many blocks? this way, one block takes about 7 min
    
    setup.ntrials           = round(setup.trialrep/setup.nblocks); % total nr of trials per block
    assert(mod(setup.ntrials, 1)==0);
    
    % set a high motion difference for the example
    setup.threshold         = .3;
    
    % timing
    setup.fixtime           = (.5 + .5*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
    setup.viewingtime       = 1; % viewing duration in seconds (fixed in this script, or maximum viewing duration in RT paradigm
    setup.nframes           = ceil(setup.viewingtime*window.frameRate); %number of frames the stimulus is displayed
    setup.intervaltime      = (.3 + .4*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
    setup.resptime          = 1; % maximum time for response
    setup.pupilreboundtime  = (.5 + .5 * rand(setup.nblocks, setup.ntrials)); %1.5-2.5s after resp
    setup.pupilreboundtime2 = (.5 + .5 * rand(setup.nblocks, setup.ntrials));% 2-2.5s after fb
    
    % direction depends on the participant
    dots.direction = 90;
    
    % decrease or increase from baseline coherence
    setup.increments         = [1 -1];
    setup.increment          = [1 -1 1 -1];
    
    % for each trial, compute the actual coherence (from the baseline, individual threshold and the sign of the increment).
    dots.coherence      = setup.baselinecoh + setup.increment.* setup.threshold;
    
    %% preload all the dot coordinates
    coord.fix           = nan(setup.nblocks, setup.ntrials, ceil(max(setup.fixtime(:))*window.frameRate), 2, dots.nDots);
    coord.ref           = nan(setup.nblocks, setup.ntrials, setup.nframes, 2, dots.nDots);
    coord.interval      = nan(setup.nblocks, setup.ntrials, ceil(max(setup.intervaltime(:))*window.frameRate), 2, dots.nDots);
    coord.stim          = nan(setup.nblocks, setup.ntrials, setup.nframes, 2, dots.nDots);
    coord.resp          = nan(setup.nblocks, setup.ntrials, ceil(max(setup.resptime(:))*window.frameRate), 2, dots.nDots);
    
    for block = 1:setup.nblocks,
        for trial = 1:setup.ntrials,
            
            % preload all the dot coordinates before starting the trial
            coord.fix(block, trial, :, :, :)        = dots_noise(dots, ceil(max(setup.fixtime(:))*window.frameRate));
            coord.ref(block, trial, :, :, :)        = dots_refstim(setup, window, dots, block, trial);
            coord.interval(block, trial, :, :, :)   = dots_noise(dots, ceil(max(setup.intervaltime(:))*window.frameRate));
            coord.stim(block, trial, :, :, :)       = dots_limitedlifetime(setup, window, dots, block, trial);
            coord.resp(block, trial, :, :, :)       = dots_noise(dots, ceil(max(setup.resptime(:))*window.frameRate));
            
        end
    end
    
    save(sprintf('%s/Data/MotionEnergy/motionStrength_timecourse.mat', mypath), 'window', 'coord', 'setup', 'dots');
    clear all; close all; clc;
    load(sprintf('%s/Data/MotionEnergy/motionStrength_timecourse.mat', mypath));
    
    %% RUN THE ACTUAL MOTION ENERGY FILTER ON THIS
    
    % to avoid window UI bugs
    display.dist        = window.dist;
    display.res         = window.res;
    display.width       = window.width;
    display.frameRate   = window.frameRate;
    display.center      = window.center;
    clear sound audio
    
    % run this once to find optimal fft algorithm
    fftw('planner', 'exhaustive');
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute some general things for this file
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfg.ppd             = deg2pix(display, 1);  % pixel per degree
    cfg.srange          = [-0.7 0.7];           % predetermine size the filter  space will have
    cfg.frameRate       = display.frameRate;
    
    % to avoid bugs with frameRates that are slightly below 60 Hz
    if cfg.frameRate < 60 && round(cfg.frameRate)==60,
        cfg.frameRate = 60.1;
    end
    
    cfg.trange          = [0 0.2];              % temporal range
    cfg.filtsize        = [ceil(diff(cfg.srange)*cfg.ppd) ...
        ceil(diff(cfg.srange)*cfg.ppd) ...
        ceil(diff(cfg.trange)*cfg.frameRate)];
    
    % pad with some pixels for valid convolution
    cfg.stimpad         = 50;
    
    setup.nframes       = 270; % for trial 1;
    
    cfg.stimsize        = [2*dots.radius+1+cfg.stimpad 2*dots.radius+1+cfg.stimpad setup.nframes]; %  predetermine the size the cloud of dots will have
    cfg.n_convolution   = cfg.stimsize + cfg.filtsize - 1;  % 'full'  size of the convolution
    cfg.nextpow2        = 2.^nextpow2(cfg.n_convolution);
    cfg.validsize       = cfg.stimsize - cfg.filtsize + 1;  % 'valid' size of convolution
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate filters
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    direc                 = dots.direction;
    counterdir            = dots.direction+90;
    if counterdir > 360, counterdir = counterdir - 360; end
    theta                 = [direc counterdir]; % only those two directions
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. CREATE SPATIAL AND TEMPORAL FILTERS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % all parameters from Kiani et al., http://www.jneurosci.org/content/28/12/3017.full
    
    % time and space axes
    x = cfg.srange(1) : 1/cfg.ppd : cfg.srange(2);  % range in space, in degree, spaced by the size in degree of one pixel
    y = cfg.srange(1) : 1/cfg.ppd : cfg.srange(2);  % range in space, in degree, spaced by the size in degree of one pixel
    t = cfg.trange(1) : 1/cfg.frameRate: cfg.trange(2);   % range in time, in s, spaced by framerate
    
    % important: use an uneven nr of points, easier for the later convolution
    assert(mod(numel(x), 2) == 1, 'x should have an uneven nr of samples');
    assert(mod(numel(y), 2) == 1, 'y should have an uneven nr of samples');
    assert(mod(numel(t), 2) == 1, 't should have an uneven nr of samples');
    
    % check if the filter has the size I thought it would get
    assert(all([length(x) length(y) length(t)] == cfg.filtsize), 'filter does not have prespecified size');
    
    % constants of the spatial filters
    sc      = 0.35;
    sg      = 0.05;
    
    % create spatial mesh
    [xmesh, ymesh] = meshgrid(x,y);
    
    % constants of the temporal filters
    k       = 60;
    nslow   = 3;
    nfast   = 5;
    
    % TEMPORAL FUNCTIONS
    % put them in the third dimension for the outer product
    g1(1,1,:) = (k*t).^nslow .* exp(-k*t) .* (1/factorial(nslow) - ((k*t).^2) / factorial(nslow + 2));
    g2(1,1,:) = (k*t).^nfast .* exp(-k*t) .* (1/factorial(nfast) - ((k*t).^2) / factorial(nfast + 2));
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make a separate filter for each direction
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for thistheta = theta,
        
        % rotate the meshgrid
        yprime  = - xmesh .* sind(thistheta) + ymesh .* cosd(thistheta);
        xprime  =   xmesh .* cosd(thistheta) + ymesh .* sind(thistheta);
        
        % SPATIAL FUNCTIONS two fourth order cauchy functions
        % transpose to match my unit circle directionality
        alpha      = atand(xprime ./ sc);
        f1rot      = (cosd(alpha).^4 .* cosd(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
        f2rot      = (cosd(alpha).^4 .* sind(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
        
        % these two linear filters are in space-time quadrature
        filt1   = bsxfun(@times, f1rot, g1) + bsxfun(@times, f2rot, g2);
        filt2   = bsxfun(@times, f2rot, g1) - bsxfun(@times, f1rot, g2);
        
        % save the filters to output
        % already run the filter fft!
        if 0,
            filters.one(:, :, :, find(thistheta==theta)) = fftn(single(filt1), cfg.n_convolution);
            filters.two(:, :, :, find(thistheta==theta)) = fftn(single(filt2), cfg.n_convolution);
        else
            filters.(['theta' num2str(find(thistheta==theta))]).one = filt1;
            filters.(['theta' num2str(find(thistheta==theta))]).two = filt2;
        end
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load stimuli and run convolution operation
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for trial = 1; %:setup.ntrials,
        
        stimulus = cat(1, ...
            squeeze(coord.fix(1, trial, :, :, :)), ...
            squeeze(coord.ref(1, trial, :, :, :)), ...
            squeeze(coord.interval(1, trial, :, :, :)), ...
            squeeze(coord.stim(1, trial, :, :, :)), ...
            squeeze(coord.resp(1, trial, :, :, :))) ;
        assert(~any(isnan(stimulus(:))));
        
        % take this from the setup struct
        xres    = display.res.width;
        yres    = display.res.height;
        nfr     = size(stimulus, 1);
        
        % preallocate
        stimrep = zeros(xres, yres, nfr);
        
        for f = 1:nfr,
            
            posx = squeeze(stimulus(f, 1, :));
            posy = squeeze(stimulus(f, 2, :));
            
            % normalize by the range of values, so that all position estimates are
            % round to the nearest pixel to create matrix
            posx = round(posx + display.center(1));
            posy = round(posy + display.center(2));
            
            % put those coordinates in the stim representation matrix
            for i = 1:length(posx), % for each dot
                % Y X T
                stimrep(posx(i), posy(i), f) = 1; % put a pixel in the matrix
            end
        end
        
        % continue with only those parts of the stimulus in
        % the x and y dir that contain the cloud of dots
        x2use = [display.center(1)-dots.radius-cfg.stimpad/2 : display.center(1)+dots.radius+cfg.stimpad/2];
        y2use = [display.center(2)-dots.radius-cfg.stimpad/2 : display.center(2)+dots.radius+cfg.stimpad/2];
        
        stimrep = stimrep(x2use, y2use, :);
        
        % get the fft of the stimulus
        stim_fft = fftn(stimrep, cfg.n_convolution);
        
        % clear stimulus stimrep
        % determine the size the output will have
        for s = 1:3,
            size2use(s,:) = [cfg.stimsize(s)-cfg.validsize(s)+1  ...
                cfg.n_convolution(s) - (cfg.stimsize(s)-cfg.validsize(s)) ];
        end
        
        %% 2. FILTER AT EACH THETA
        for thistheta = theta,
            
            % run multiplication in the frequency domain
            % this is where the bulk of the computation happens
            
            %% speed test
            
            cfg.n_convolution = size(stimrep);
            resp1       = ifftn(fftn(stimrep, cfg.n_convolution) .* fftn(filters.(['theta' num2str(find(thistheta==theta))]).one, cfg.n_convolution), cfg.n_convolution);
            resp2       = ifftn(stim_fft .* filters.two(:, :, :, find(thistheta==theta)), cfg.n_convolution);
            
            % use only valid part of the result, slightly smaller than the size of the input
            resp1       = resp1(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
            resp2       = resp2(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
            
            % sum and square the results of the two filters in quadrature, see Adelson & Bergen
            energy      = (resp1.^2 + resp2.^2);
            
            % collapse over the x and y directions, just give the time output
            % take the square root, to normalize the responses
            motionenergy(find(thistheta == theta), :) = single(sqrt(squeeze(sum(sum(energy)))))';
            
        end % theta
    end
    
    motionstrength = squeeze(motionenergy(1, :) - motionenergy(2,:));
    save(sprintf('%s/Data/MotionEnergy/motionTimecourse.mat', mypath), 'motionstrength');
    
else
    % get the file we already have
    load(sprintf('%s/Data/MotionEnergy/motionTimecourse.mat', mypath));
    window.frameDur = 1/60;
    
    timeaxis = 0:window.frameDur:length(motionstrength)*window.frameDur;
    plot(timeaxis(1:end-1),motionstrength, 'k');
    axis tight; xlabel('Time (s)'); ylabel('Motion strength'); box off;
    ylim([-1 11]);
    
end

end

function pix = deg2pix(display, ang)

% Converts visual angles in degrees to pixels.
%
% Inputs:
% display.dist (distance from screen (cm))
% display.width (width of screen (cm))
% display.res (number of pixels of display in horizontal direction)
% ang (visual angle)
% Warning: assumes isotropic (square) pixels
%
% Written 11/1/07 gmb zre
% Adapted by Anne Urai, August 2013

pixSize = display.width/display.res.width;   %cm/pix

sz = 2*display.dist*tan(pi*ang/(2*180));  %cm

pix = round(sz/pixSize);   %pix

return

end

function ang = pix2deg(window,pix)
%angle = pix2angle(window,pix)
%
%converts monitor pixels into degrees of visual angle.
%
%Inputs:
%window.dist (distance from screen (cm))
%window.width (width of screen (cm))
%window.resolution (number of pixels of window in horizontal direction)
%
%ang (visual angle)
%
%Warning: assumes isotropic (square) pixels

%Written 11/1/07 gmb zre

%Calculate pixel size
pixSize = window.width/window.res.width;   %cm/pix

sz = pix*pixSize;  %cm (duh)

ang = 2*180*atan(sz/(2*window(1).dist))/pi;


return
end


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

function [ stimuli] = dots_limitedlifetime(setup, window, dots, block, trial)

% Generates a dot patch according to a limited lifetime
% algorithm, see:
% Pilly, P. K., & Seitz, A. R. (2009).
% What a difference a parameter makes: a psychophysical comparison
% of random dot motion algorithms. Vision Research, 49(13), 1599?612.

% IMPORTANT direction of motion is clockwise from vertical, so not the unit
% circle!
%--------------------------------------------------------------------------

NVAR        = dots.nvar; % interleaved sequence
LIFETIME    = dots.lifetime;
NDOTS       = dots.nDots; % 774, 3097

if length(dots.direction) > 1,
    DIRECTION   = dots.direction(block, trial);
else
    DIRECTION   = dots.direction;
end

COH         = dots.coherence(block,trial);
SPEED       = dots.speed; % downward, as a function of speed (in degrees per second) and framerate)
NFR         = ceil(setup.nframes/NVAR);
RADIUS      = dots.radius;
INNER       = dots.innerspace; %inner circle doesn't contain dots

for var = 1:NVAR, %create separate variants
    
    % generate random starting points within a circular aperture
    rad     = RADIUS*sqrt(rand(NDOTS,1));
    theta   = 2*pi*rand(NDOTS,1);
    pos     = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates
    
    % assign each dot a scalar that indicates how long it has been
    % 'alive' (that is, a signal dot)
    % Each dot will have a integer value 'life' which is how many frames the
    % dot has been going.  The starting 'life' of each dot will be a random
    % number between 0 and dots.lifetime-1 so that they don't all 'die' on the
    % same frame:
    
    life    = ceil(rand(1,NDOTS)*LIFETIME);
    
    for frameNum = 1:NFR,
        
        % define the noise dots
        noisedots           = round(NDOTS*(1-COH)); %nr of noisedots based on the coherence
        t1                  = rand(NDOTS,1);
        [t1,t2]             = sort(t1);
        noiseindex          = t2(1:noisedots); %random subset of dots
        % define signal dots
        signalindex         = t2(noisedots+1:end); %the dots that are signal
        
        % move signal dots with a certain speed in the right direction
        pos(signalindex,:)  = [SPEED.*cos(DIRECTION*pi/180)/window.frameRate+pos(signalindex,1), ...
            SPEED.*sin(DIRECTION*pi/180)/window.frameRate+pos(signalindex,2)]; %convert to cartesian coordinates
        
        % replot the noisedots somewhere in the aperture
        rad                 = RADIUS*sqrt(rand(noisedots,1));
        theta               = 2*pi*rand(noisedots,1);
        [pos(noiseindex, 1), pos(noiseindex,2)] = pol2cart(theta, rad);
        
        %increment the 'life' of each dot
        life                = life+1;
        
        %find the 'dead' dots
        deadindex           = mod(life,LIFETIME)==0;
        deaddots            = length(find(deadindex==1));
        
        %replace the positions of the dead dots to a random location
        rad                 = RADIUS*sqrt(rand(deaddots,1));
        theta               = 2*pi*rand(deaddots,1);
        [pos(deadindex, 1), pos(deadindex,2)] = pol2cart(theta, rad);
        
        %find the dots that have left the aperture
        outindex            = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) >= RADIUS); %index dots that are outside the aperture
        
        % wrap them around in the direction where they came from
        [theta, rad]        = cart2pol(pos(outindex, 1), pos(outindex,2));
        theta               = theta + pi; %move to the other side of the circle
        rad                 = RADIUS*ones(length(rad),1);
        [pos(outindex, 1), pos(outindex,2)] = pol2cart(theta, rad);
        
        %find the dots that are too close to the fixation
        innerindex          = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) <= INNER);
        [theta, rad]        = cart2pol(pos(innerindex, 1), pos(innerindex, 2));
        rad                 = INNER + (RADIUS - INNER)*sqrt(rand(length(innerindex),1)); %random radius
        theta               = theta + pi;
        [pos(innerindex, 1), pos(innerindex,2)] = pol2cart(theta, rad);
        
        % save the positions per variant
        temp(var,frameNum, :, :)  = pos';
    end
end

stimuli = squeeze(reshape(temp, 1, NFR*NVAR, 2, NDOTS));
stimuli = stimuli(1:setup.nframes, :, :);

end

function [ stimuli] = dots_refstim(setup, window, dots, block, trial, stimuli)

% Generates a dot patch according to a limited lifetime
% algorithm, see:
% Pilly, P. K., & Seitz, A. R. (2009).
% What a difference a parameter makes: a psychophysical comparison
% of random dot motion algorithms. Vision Research, 49(13), 1599?612.
%--------------------------------------------------------------------------

NVAR        = dots.nvar; % interleaved sequence
LIFETIME    = dots.lifetime;
NDOTS       = dots.nDots; % 774, 3097
DIRECTION   = dots.direction;
SPEED       = dots.speed; % downward, as a function of speed (in degrees per second) and framerate)
NFR         = ceil(setup.nframes/NVAR);
RADIUS      = dots.radius;
INNER       = dots.innerspace; %inner circle doesn't contain dots
COH         = setup.baselinecoh;

for var = 1:NVAR, %create separate variants
    
    % generate random starting points within a circular aperture
    rad     = RADIUS*sqrt(rand(NDOTS,1));
    theta   = 2*pi*rand(NDOTS,1);
    pos     = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates
    
    % assign each dot a scalar that indicates how long it has been
    % 'alive' (that is, a signal dot)
    % Each dot will have a integer value 'life' which is how many frames the
    % dot has been going.  The starting 'life' of each dot will be a random
    % number between 0 and dots.lifetime-1 so that they don't all 'die' on the
    % same frame:
    
    life    = ceil(rand(1,NDOTS)*LIFETIME);
    
    for frameNum = 1:NFR,
        
        % define the noise dots
        noisedots           = round(NDOTS*(1-COH)); %nr of noisedots based on the coherence
        t1                  = rand(NDOTS,1);
        [t1,t2]             = sort(t1);
        noiseindex          = t2(1:noisedots); %random subset of dots
        % define signal dots
        signalindex         = t2(noisedots+1:end); %the dots that are signal
        
        % move signal dots with a certain speed in the right direction
        pos(signalindex,:)  = [SPEED.*cos(DIRECTION*pi/180)/window.frameRate+pos(signalindex,1), ...
            SPEED.*sin(DIRECTION*pi/180)/window.frameRate+pos(signalindex,2)]; %convert to cartesian coordinates
        
        % replot the noisedots somewhere in the aperture
        rad                 = RADIUS*sqrt(rand(noisedots,1));
        theta               = 2*pi*rand(noisedots,1);
        pos(noiseindex,:)   = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates
        
        %increment the 'life' of each dot
        life                = life+1;
        
        %find the 'dead' dots
        deadindex           = mod(life,LIFETIME)==0;
        deaddots            = length(find(deadindex==1));
        
        %replace the positions of the dead dots to a random location
        rad                 = RADIUS*sqrt(rand(deaddots,1));
        theta               = 2*pi*rand(deaddots,1);
        pos(deadindex,:)    = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates
        
        %find the dots that have left the aperture
        outindex            = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) >= RADIUS); %index dots that are outside the aperture
        % wrap them around in the direction where they came from
        [theta, rad]        = cart2pol(pos(outindex, 1), pos(outindex,2));
        theta               = theta + pi; %move to the other side of the circle
        rad                 = RADIUS*ones(length(rad),1);
        pos(outindex, :)    = [rad.*cos(theta) rad.*sin(theta)]; %move back to the top, only change y coordinate
        
        %find the dots that are too close to the fixation
        innerindex          = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) <= INNER);
        [theta, rad]        = cart2pol(pos(innerindex, 1), pos(innerindex, 2));
        rad                 = INNER + (RADIUS - INNER)*sqrt(rand(length(innerindex),1)); %random radius
        theta               = theta + pi;
        pos(innerindex, :)  = [rad.*cos(theta) rad.*sin(theta)]; %move back to the top, only change y coordinate
        
        % save the positions per variant
        temp(var,frameNum, :, :)  = pos';
    end
end

stimuli = squeeze(reshape(temp, 1, NFR*NVAR, 2, NDOTS));
stimuli = stimuli(1:setup.nframes, :, :);

end
