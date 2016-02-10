function [] =  a2_MotionEnergy(sj)
% for (almost) all subjects, the coordinates of each trial's dot were saved
% in the behavioural files. Using the method described by Adelson & Bergen
% and Kiani et al, these stimulus movies are filtered to compute
% single-trial timecourses of fluctuations in motion energy
%
% Anne Urai, 2015

if ischar(sj), sj = str2double(sj); end

global mypath;
mypath = '~/Data/pupilUncertainty';

% create logfile (handy when running on the cluster, the script will find
% which subject to work on by itself)
cd(sprintf('%s/Data/P%02d/', mypath, sj));
system(sprintf('touch P%02d_motionEstarted.log', sj));

clear sessions;
% check which sessions to use
s = dir('S*');
s = {s(:).name};
for i = 1:length(s), sessions(i) = str2num(s{i}(2)); end

for session = sessions,
    
    % preallocate the motion energy output
    mdat.int1     = single(nan(10, 50));
    mdat.int2     = single(nan(10, 50));
    mdat.strength = single(nan(10, 50));
    mdat.stim     = single(nan(10, 50));
    mdat.response = single(nan(10, 50));
    mdat.correct  = single(nan(10, 50));
    
    trllength       = 18;
    mdat.timecourse = single(nan(10,50,2,2,trllength)); % block trl int theta time
    
    % ==================================================================
    % some subjects didnt do all blocks, manually correct
    % ==================================================================
    
    if sj == 2 && session == 1,
        blocks = 4:9;
    elseif sj == 4 && session == 1,
        blocks = 6:10; % first 5 blocks, lost dot coords
    elseif sj == 17 && session == 1,
        blocks = 1:5;
    elseif sj == 15 && session == 3,
        blocks = 1:5;
    elseif sj == 15 && session == 6,
        blocks = 1:5;
    elseif sj == 12 && session == 1,
        blocks = 7:10;
    elseif sj == 8 && session == 4,
        blocks = 6:10;
    elseif sj == 5 && session == 1,
        blocks = [8 9 10 4:10];
    else
        blocks = 1:10;
    end
    
    % some weirdness in the blockcount, correct for this
    if sj == 5 && session == 1,
        blockidx = 1:10;
    else
        blockidx = blocks;
    end
    
    blockcnt = 0;
    for b = 1:length(blocks),
        iblock = blocks(b);
        blockcnt = blockcnt + 1;
        
        % ==================================================================
        % get the files we need
        % ==================================================================
        
        clearvars -except subjects sj sessions session blocks iblock blockidx blockcnt mdat mypath
        begin = tic;
        
        if sj > 14,
            % for all the SJs OJay measured, there are separate coord files
            files = dir(sprintf('%s/Data/P%02d/Behav/Dots_P%d_s%d_b%d_20*.mat', mypath, sj, sj, session, iblock));
            assert(length(files)==1);
            load(sprintf('%s/Data/P%02d/Behav/%s', mypath, sj, files.name));
            fprintf('%s/Data/P%02d/Behav/%s \n',  mypath, files.name);
            
        else
            % get the one behav file for this session, includes the coords
            files = dir(sprintf('%s/Data/P%02d/Behav/P%d_s%d_20*.mat', mypath, sj, sj, session));
            if length(files) == 1,
                load(sprintf('%s/Data/P%02d/Behav/%s', mypath, sj, files.name));
                fprintf('%s/Data/P%02d/Behav/%s \n', mypath, sj, files.name);
                
            else % if there are multiple files
                for f = 1:length(files),
                    load(sprintf('%s/Data/P%02d/Behav/%s',  mypath, sj, files(f).name));
                    fprintf('%s/Data/P%02d/Behav/%s \n',  mypath, sj, files(f).name);
                    % check if there are responses in this block
                    try
                        if ~isnan(nanmean(results.response(iblock, :))) && ...
                                length(find(~isnan(results.response(iblock, :)))) > 20,
                            break % only continue if this block has all the responses we need
                        end
                    end
                end % try files
            end % try files
        end % sj nr
        
        % something weird with SJ 5...
        if sj == 5 && session == 1 && blockcnt == 3,
            cd(sprintf('%s/Data/P05/Behav', mypath));
            system('mv P5_s1_2014-02-24_15-07-33.mat first3_P5_s1_2014-02-24_15-07-33.mat');
        end
        
        fprintf('starting P%02d_s%d_b%d \n', sj, session, iblock);
        
        % ==================================================================
        % change the window structure into display, window is a
        % matlab function
        % ==================================================================
        
        % to avoid window UI bugs
        display.dist        = window.dist;
        display.res         = window.res;
        display.width       = window.width;
        display.frameRate   = window.frameRate;
        display.center      = window.center;
        clear sound audio
        
        % run this once to find optimal fft algorithm
        fftw('planner', 'exhaustive');
        
        % ==================================================================
        % compute features of the filters based on this file
        % ==================================================================
        
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
        cfg.stimsize        = [2*dots.radius+1+cfg.stimpad 2*dots.radius+1+cfg.stimpad setup.nframes]; %  predetermine the size the cloud of dots will have
        cfg.n_convolution   = cfg.stimsize + cfg.filtsize - 1;  % 'full'  size of the convolution
        cfg.nextpow2        = 2.^nextpow2(cfg.n_convolution);
        cfg.validsize       = cfg.stimsize - cfg.filtsize + 1;  % 'valid' size of convolution
        
        % ==================================================================
        % find directions in which to filter
        % ==================================================================
        
        direc                 = dots.direction;
        counterdir            = dots.direction + 180;
        if counterdir > 360, counterdir = counterdir - 360; end
        
        theta                 = [direc counterdir]; % only those two directions
        
        % ==================================================================
        % 1. CREATE SPATIAL AND TEMPORAL FILTERS
        % ==================================================================
        
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
        
        % ==================================================================
        % make a separate filter for each direction
        % ==================================================================
        
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
            filters.one(:, :, :, find(thistheta==theta)) = fftn(single(filt1), cfg.n_convolution);
            filters.two(:, :, :, find(thistheta==theta)) = fftn(single(filt2), cfg.n_convolution);
            
            if 0,
                
                % ==================================================================
                % can plot an overview of what the filters look like
                % ==================================================================
                
                figure; colormap bone;
                subplot(441); imagesc(x, y, f1rot); xlabel('x'); ylabel('y'); title('f1');
                subplot(445); imagesc(x, y, f2rot); xlabel('x'); ylabel('y'); title('f2');
                
                subplot(4,4,[2 6]); plot(t, squeeze(g1), t, squeeze(g2)); title('g1 g2'); xlabel('time');
                
                % freq spectrum of filters
                subplot(4,4,3); imagesc(x, t, squeeze(sum((filter.one),2)) .* squeeze(sum((filter.two),2)));
                %xlabel('time'); ylabel('x');
                title('freq along xdir');
                subplot(4,4,7); imagesc(y, t, squeeze(sum((filter.one),1)) .* squeeze(sum((filter.two),1)));
                %xlabel('time'); ylabel('y');
                title('freq along ydir');
                
                % collapsed over x dir
                subplot(4,4,9); imagesc(y, t, squeeze(sum(filt1,1)));
                title('filt1'); xlabel('y'); ylabel('t');
                subplot(4,4,13); imagesc(y, t, squeeze(sum(filt2,1)));
                title('filt2'); xlabel('y'); ylabel('t');
                
                % collapse over y dir
                subplot(4,4,10); imagesc(x, t, squeeze(sum(filt1, 2)));
                title('filt1'); xlabel('x'); ylabel('t');
                subplot(4,4, 14); imagesc(x, t, squeeze(sum(filt2,2)));
                title('filt2'); xlabel('x'); ylabel('t');
                
                % collapse over t dir
                subplot(4,4,11); imagesc(x, y, squeeze(sum(filt1, 3)));
                title('filt1'); xlabel('x'); ylabel('y');
                subplot(4,4, 15); imagesc(x, y, squeeze(sum(filt2,3)));
                title('filt2'); xlabel('x'); ylabel('y');
                
                suplabel([ 'Rotation ' num2str(thistheta) ' degrees'], 't');
                saveas(gcf, sprintf('~/Data/MotionEnergy/filters/filterimg_%d.eps', thistheta), 'epsc');
                close all;
            end
        end
        
        fprintf('preparing filters took %.2f seconds \n', toc(begin));
        
        % ==================================================================
        % convert stim coordinates into movie matrix
        % ==================================================================
        
        for trial = 1:setup.ntrials,
            motionenergy = nan(2, length(theta), cfg.validsize(3));
            
            for int = 1:2, % loop over reference and test interval
                
                % ==================================================================
                % get the coordinates that we need, this is quite messy because the format
                % is different for different subjects
                % ==================================================================
                
                % for those subjects where there is no refstim
                % coords saved, continue
                if sj < 15 && ~exist('refstim', 'var') && int == 1,
                    continue;
                end
                
                % stimulus should be of dimensions: nframes, [x y], ndots
                % why on earth did I saved this in so many different ways
                if sj > 14,
                    % subjects that OJay measured, coord structure
                    switch int
                        case 1,
                            stimulus = squeeze(coord.ref(1, trial, :, :, :));
                        case 2,
                            stimulus = squeeze(coord.stim(1, trial, :, :, :));
                    end
                elseif iscell(stim),
                    % for some early participants, yet another format
                    switch int
                        case 1,
                            error('this is not supposed to happen');
                        case 2,
                            for f = 1:size(stim{iblock, trial}, 3),
                                stimulus(f, 1, :) = squeeze(stim{iblock, trial, :}(iblock, trial, f).x);
                                stimulus(f, 2, :) = squeeze(stim{iblock, trial, :}(iblock, trial, f).y);
                            end
                    end
                else
                    % sj < 15, 5d matrix in each behav session file
                    switch int
                        case 1,
                            stimulus = squeeze(refstim(iblock, trial, :, :, :));
                        case 2,
                            stimulus = squeeze(stim(iblock, trial, :, :, :));
                    end
                end
                
                assert(~any(isnan(stimulus(:))));
                
                % take this from the setup struct
                xres    = display.res.width;
                yres    = display.res.height;
                nfr     = setup.nframes;
                
                % ==================================================================
                % generate a frame-by-frame image movie of the stim
                % ==================================================================
                
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
                assert(all(size(stimrep) == cfg.stimsize), 'stimulus size not correctly computed');
                
                % get the fft of the stimulus
                stim_fft = fftn(stimrep, cfg.n_convolution);
                
                clear stimulus stimrep
                % determine the size the output will have
                for s = 1:3,
                    size2use(s,:) = [cfg.stimsize(s)-cfg.validsize(s)+1  ...
                        cfg.n_convolution(s) - (cfg.stimsize(s)-cfg.validsize(s)) ];
                end
                
                % ==================================================================
                % RUN FILTER AT EACH THETA
                % ==================================================================
                
                for thistheta = theta,
                    
                    % run multiplication in the frequency domain
                    % this is where the bulk of the computation happens
                    
                    resp1       = ifftn(stim_fft .* filters.one(:, :, :, find(thistheta==theta)), cfg.n_convolution);
                    resp2       = ifftn(stim_fft .* filters.two(:, :, :, find(thistheta==theta)), cfg.n_convolution);
                    
                    % use only valid part of the result, slightly smaller than the size of the input
                    resp1       = resp1(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
                    resp2       = resp2(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
                    
                    % sum and square the results of the two filters in quadrature, see Adelson & Bergen
                    energy      = (resp1.^2 + resp2.^2);
                    
                    % collapse over the x and y directions, just give the time output
                    % take the square root, to normalize the responses
                    motionenergy(int, find(thistheta == theta), :) = single(sqrt(squeeze(sum(sum(energy)))))';
                    
                end % theta
            end %int
            
            % ==================================================================
            % TRANSFORM FILTER TIMECOURSES INTO SINGLE-TRIAL SCALARS
            % ==================================================================
            
            if sj == 5 && session == 1,
                % take the difference between the two intervals for this trial
                motiondiffRef                 = squeeze(motionenergy(1, 1, :) - motionenergy(1, 2, :));
                mdat.int1(blockcnt, trial)    = mean(motiondiffRef);
                if isnan(mdat.int1(blockcnt, trial)),
                    % if there were no coordinates, use the mean...
                    mdat.int1(blockcnt, trial)  = 12.6;
                end
                
                % every dataset should have this, even if there is no refstim
                motiondiffStim = squeeze(motionenergy(2, 1, :) - motionenergy(2, 2, :));
                mdat.int2(blockcnt, trial)      = mean(motiondiffStim);
                
                % compute a scalar motionstrength measure
                mdat.strength(blockcnt, trial)  = mdat.int2(blockcnt, trial) - mdat.int1(blockcnt, trial) ;
                
                % save the full timecourse as well
                mdat.timecourse(blockcnt, trial, :, :, 1:cfg.validsize(3)) = motionenergy;
                
            else
                
                % take the difference between the two intervals for this trial
                motiondiffRef                 = squeeze(motionenergy(1, 1, :) - motionenergy(1, 2, :));
                mdat.int1(iblock, trial)      = mean(motiondiffRef);
                if isnan(mdat.int1(iblock, trial)),
                    % if there were no coordinates, use the mean...
                    mdat.int1(iblock, trial)  = 12.6;
                end
                
                % every dataset should have this, even if there is no refstim
                motiondiffStim = squeeze(motionenergy(2, 1, :) - motionenergy(2, 2, :));
                mdat.int2(iblock, trial)      = mean(motiondiffStim);
                
                % compute a scalar motionstrength measure
                mdat.strength(iblock, trial)  = mdat.int2(iblock, trial) - mdat.int1(iblock, trial) ;
                
                % save the full timecourse as well
                mdat.timecourse(iblock, trial, :, :, 1:cfg.validsize(3)) = motionenergy;
                
            end
        end % trial
        
        % ==================================================================
        % save some additional behavioural data to match this to pupil-derived
        % datafiles later
        % ==================================================================
        
        if sj == 15 && session == 6,
            % add extra info that will make it easier to match with pupil files
            mdat.stim(blockidx(iblock), :)         = dots.coherence(iblock+5, :) - 0.7;
            mdat.response(blockidx(iblock), :)     = results.response(iblock+5, :);
            mdat.correct(blockidx(iblock), :)      = results.correct(iblock+5, :);
        elseif sj == 5 && session == 1,
            disp('using blockcnt');
            % add extra info that will make it easier to match with pupil files
            mdat.stim(blockcnt, :)         = dots.coherence(iblock, :) - 0.7;
            mdat.response(blockcnt, :)     = results.response(iblock, :);
            mdat.correct(blockcnt, :)      = results.correct(iblock, :);
        else
            % add extra info that will make it easier to match with pupil files
            mdat.stim(iblock, :)         = dots.coherence(iblock, :) - 0.7;
            mdat.response(iblock, :)     = results.response(iblock, :);
            mdat.correct(iblock, :)      = results.correct(iblock, :);
        end
    end % block
    
    save(sprintf('%s/Data/MotionEnergy/motionenergy_P%02d_s%d.mat', mypath, sj, session), '-mat', 'mdat');
end

system(sprintf('touch P%02d_motionEfinished.log', sj));

end% function end

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

