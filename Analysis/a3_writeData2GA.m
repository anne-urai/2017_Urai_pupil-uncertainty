function [] = a3_writeData2GA()
% write one grand average of pupil values for all subjects, on which we can
% do the regression timecourse and define the window we will work with.
%
% Anne Urai, 2015

global mypath;

subjects = 1:27;
for sj = fliplr(subjects),
    tic;
    
    clearvars -except sj subjects alldat pupilgrandavg;
    % choose between 2 and 3
    load(sprintf('%s/Data/P%02d_alleye.mat', mypath, sj));
    
    % check which sessions to use
    if sj == 15, sessions = 1:7;
    else sessions = 1:6;
    end
    
    % ==================================================================
    % GET ALL THE FILTERED MOTIONENERGY FOR THIS SUBJECT
    % ==================================================================
    
    clear mdats mdat
    for session = sessions,
        load(sprintf('%s/Data/MotionEnergy/motionenergy_P%02d_s%d.mat', mypath, sj, session));
        
        % transform into table
        mdat = structfun(@transpose, mdat, 'uniformoutput', 0);
        mdat = structfun(@(x) reshape(x, 500, 1), mdat, 'uniformoutput', 0);
        mdat = struct2table(mdat);
        
        % add info about trialnum, blocknum, sessionnum
        mdat.trialnum   = repmat(1:50', 1, height(mdat)/50)';
        mdat.blocknum   = reshape(repmat(1:10, 50, 1), height(mdat), 1);
        mdat.session    = session*ones(size(mdat.int1));
        mdats{session}  = mdat;
    end
    
    % merge all 6 sessions
    mdats          = cat(1, mdats{:});
    trl            = data.trialinfo;
    trl(:, 15)     = zeros(size(trl(:,14)));
    
    % correct for a weird mistake
    if sj == 3,
        trl(find(abs(trl(:, 4) - 0.01) < 0.00001), 4) = 0.0125;
    end
    
    % ==================================================================
    % match the single-trial motionenergy to the trl matrix
    % ==================================================================
    
    for t = 1:length(trl),
        
        % find the matching idx
        thist = find(trl(t,12) == mdats.trialnum ...
            & trl(t,13) == mdats.blocknum & trl(t,14) == mdats.session);
        
        try
            % check that this makes sense
            assert(~isempty(thist));
            assert(trl(t,3) .* trl(t,4) - mdats.stim(thist) < 0.001);
            % then add to the trlinfo
            trl(t,15) = mdats.strength(thist);
        catch
            disp('could not match trials!');
            trl(t, 15) = 0; % fill with nan for now
        end
    end
    
    % for some subjects, the dot coordinates were lost
    if sj == 12 || sj == 4 || sj == 3 || sj == 8,
        
        % replace those missing values with the means from other trials
        wrongtrials = find(trl(:, 15) == 0);
        if sj == 3,
            wrongtrials = [find(trl(:, 14) == 1); find(trl(:, 14) == 2 & trl(:, 13) > 1)];
        end
        trl(wrongtrials, 15) = 0;
        
        % this seems to go wrong!
        fprintf('sj %d, %d trials without match found \n', sj, length(find(trl(:, 15)==0)));
        thisstim = trl(wrongtrials, 3) .* trl(wrongtrials, 4);
        stimlevels = unique(thisstim);
        
        % replace the missing values with the mean of previously filtered dots
        alldat2 = cat(1, alldat{:});
        allstim = alldat2(:, 1) .* alldat2(:, 2);
        
        meanEnergy = nan(size(stimlevels));
        for s = 1:length(stimlevels),
            meanEnergy(s) = nanmean(alldat2(find(allstim == stimlevels(s)), 4));
            trls = find(thisstim == stimlevels(s) & trl(wrongtrials, 15) == 0);
            assert(~isempty(trls));
            trl(wrongtrials(trls), 15) = meanEnergy(s);
        end
    end
    
    % check that this all worked
    assert(~any(trl(:,15)==0), 'merging motionenergy failed');
    assert(~any(isnan(trl(:,15))), 'matching motionenergy failed');
    
    fprintf('\n\nout of %d trials, %d trials not matched \n\n', length(trl(:, 15)), ...
        length(find(isnan(trl(:, 15)))));
    
    % add in an extra 'correct' field
    newcorrect = zeros(size(trl(:, 15)));
    for t = 1:length(trl),
        if sign(trl(t, 15)) == sign(trl(t, 7)),
            newcorrect(t) = 1;
        end
    end
    
    % put back into fieldtrip struct
    data.trialinfo = [trl newcorrect];
    
    % ==================================================================
    % REMOVE ANY NO-RESP TRIALS
    % ==================================================================
    
    cfg         = [];
    cfg.trials  = find(~isnan(data.trialinfo(:,7)));
    data        = ft_selectdata(cfg, data);
    
    % ==================================================================
    % use subfunction to get all the pupil info we're interested in
    % ==================================================================
    
    data.fsample          = 100; % make sure we use the resampled frequency from the pupilAnalysis
    [timelock, trialinfo] = s2b_GetIndividualData(data, sj, 0);
    
    % trialinfo matrix as it is
    newtrl         = data.trialinfo;
    RT             = (newtrl(:, 9) - newtrl(:,6)) / data.fsample;
    RT             = RT - 0.5; % the stimulus presentation was 500 ms, makes more sense to compute RT from its offset
    
    % remove sample idx
    newtrl(:, [1 2 6 9 10 11]) = [];
    
    % add in reaction times and newcorrect
    newtrl = [newtrl(:, [1:3]) newtrl(:, 9) newtrl(:, 4) RT newtrl(:, 5)...
        newtrl(:, 10) newtrl(:, [6:8]) trialinfo(:, 1)];
    
    % replace
    timelock(4).lock.trialinfo = newtrl;
    
    % also save the mat file
    for i = 1:4, timelock(i).lock = rmfield(timelock(i).lock, 'cfg'); end
    
    pupilgrandavg.timelock{find(sj==subjects)} = timelock;
    toc;
    
    % save for motion coordinate filling in
    alldat{find(sj==subjects)} = newtrl;
    
end

% ==================================================================
% write a grand average file with all the timelocked data
% ==================================================================

disp('saving timelock...');
cd mypath/Data;
mkdir GrandAverage;
savefast('%s/Data/GrandAverage/pupilgrandaverage.mat', 'pupilgrandavg');

end

function [timelock, trialinfo] = s2b_GetIndividualData(data, sj, plotme)
% in this function, the actual single-trial pupil values will be computed

pupilchan       = find(strcmp(data.label, 'EyePupil')==1);

if plotme,
    clf;
    cnt = 1;
    for session = 1:6,
        for block = 1:10,
            subplot(10, 6, cnt);
            try
                trls = find(data.trialinfo(:, 14) == session & data.trialinfo(:, 13)==block);
                trldat = cat(2, data.trial{trls});
                
                plot(trldat(pupilchan, :));
            end
            axis tight; axis off;
            cnt = cnt + 1;
        end
    end
end

% ==================================================================
% PRE-STIMULUS BASELINE
% ==================================================================

disp('baseline correcting...');
bl = nan(length(data.trial), 1);
data_blcorr = data;

for t = 1:length(data_blcorr.trial),
    
    startofref = data_blcorr.trialinfo(t,2) - data_blcorr.trialinfo(t,1);
    % use the 500ms before this (50 samples, resampled at 100Hz)
    try
        bl(t) = mean(data_blcorr.trial{t}(pupilchan, startofref-0.5*data_blcorr.fsample:startofref));
    catch
        % in case there are not enough samples before that
        bl(t) = mean(data_blcorr.trial{t}(pupilchan, 1:startofref));
    end
    data_blcorr.trial{t}(pupilchan, :) = data_blcorr.trial{t}(pupilchan, :) - bl(t); % subtract from whole trial
end

% add the baseline
trialinfo(:, 1) = bl;

cohs            = 1:5; % split by difficulty rather than physical coherence
% correct one weird mistake in sj 3
if sj == 3,
    data_blcorr.trialinfo(find(abs(data_blcorr.trialinfo(:,4)-0.01)< 0.000001), 4) = 0.0125;
end

% merge level 4 and 5
if length(cohs) == 4,
    data_blcorr.trialinfo(find(data_blcorr.trialinfo(:,5)==5), 5) = 4; % change 5 to 3
end

% ==================================================================
% TIMELOCK ALL THE 4 EPOCHS
% ==================================================================

whichlock   = {'ref', 'stim', 'resp', 'fb'};

for locking = 1:4,
    
    % lock to this offset
    switch whichlock{locking}
        case 'ref'
            offset      = data_blcorr.trialinfo(:,2) - data_blcorr.trialinfo(:,1);
            prestim     = 0.2;
            poststim    = 0.6;
        case 'stim'
            offset      = data_blcorr.trialinfo(:,6) - data_blcorr.trialinfo(:,1);
            prestim     = 0.3;
            poststim    = 0.6;
        case 'resp'
            offset      = data_blcorr.trialinfo(:,9) - data_blcorr.trialinfo(:,1);
            prestim     = 0.3;
            poststim    = 1.2;
        case 'fb'
            offset      = data_blcorr.trialinfo(:,11) - data_blcorr.trialinfo(:,1);
            prestim     = 1;
            poststim    = 2;
    end
    
    % shift the offset and do timelocking
    timelock(locking).lock = shiftoffset_timelock(data_blcorr, [], offset, prestim, poststim, data_blcorr.fsample, 0);
end

pupilchan       = find(strcmp(data.label, 'EyePupil')==1);

% ==================================================================
% DEFINE A SINGLE-TRIAL PUPIL SCALAR
% ==================================================================

% find where all of them are significant
trialinfo(:, 2)      = squeeze(nanmean(timelock(4).lock.trial(:, pupilchan, ...
    find(timelock(4).lock.time < 0 & timelock(4).lock.time > -0.250 ) ), 3));

% ==================================================================
% plot what this looks like
% ==================================================================

if plotme,
    clf;
    
    alltimelock = cat(2, squeeze(timelock(1).lock.trial(:, pupilchan, :)), ...
        squeeze(timelock(2).lock.trial(:, pupilchan, :)), ...
        squeeze(timelock(3).lock.trial(:, pupilchan, :)), ...
        squeeze(timelock(4).lock.trial(:, pupilchan, :)));
    
    % plot mean
    subplot(221);
    ph = boundedline(1:size(alltimelock, 2), squeeze(nanmean(alltimelock)), ...
        squeeze(nanstd(alltimelock)));
    
    % split by error/correct and difficulty
    thistable.correct = timelock(1).lock.trialinfo(:, 8);
    thistable.motionstrength = timelock(1).lock.trialinfo(:, 15);
    cnt = 0;
    cors = [0 1];
    for c = 1:2,
        trls = find(thistable.correct == cors(c));
        motionstrengthquantiles = quantile(abs(thistable.motionstrength(trls)), 2);
        
        for d = 1:3, % divide into 3 bins of absolute motionstrength
            
            % find those trials
            switch d
                case 1
                    trls = find(thistable.correct == cors(c) & abs(thistable.motionstrength) < motionstrengthquantiles(1));
                case 2
                    trls = find(thistable.correct == cors(c) & ...
                        abs(thistable.motionstrength) < motionstrengthquantiles(2) ...
                        &    abs(thistable.motionstrength) > motionstrengthquantiles(1) );
                case 3
                    trls = find(thistable.correct == cors(c) & ...
                        abs(thistable.motionstrength) > motionstrengthquantiles(2) );
            end
            
            cnt = cnt + 1;
            % get all timelock
            fulltimelock.mn(cnt, :) = cat(1, squeeze(nanmean(timelock(1).lock.trial(trls, pupilchan, :))), ...
                squeeze(nanmean(timelock(2).lock.trial(trls, pupilchan, :))), ...
                squeeze(nanmean(timelock(3).lock.trial(trls, pupilchan, :))), ...
                squeeze(nanmean(timelock(4).lock.trial(trls, pupilchan, :))));
            
            fulltimelock.std(cnt, :) = cat(1, squeeze(nanstd(timelock(1).lock.trial(trls, pupilchan, :))), ...
                squeeze(nanstd(timelock(2).lock.trial(trls, pupilchan, :))), ...
                squeeze(nanstd(timelock(3).lock.trial(trls, pupilchan, :))), ...
                squeeze(nanstd(timelock(4).lock.trial(trls, pupilchan, :))));
        end
    end
    
    subplot(222);
    cols = cbrewer('div', 'RdYlGn', 10);
    cols = cols([1:3 end-2:end], :);
    
    ph = boundedline(1:size(fulltimelock.mn, 2), fulltimelock.mn', ...
        permute(fulltimelock.std, [2 3 1]) ./ 4, 'cmap', cols, 'alpha');
    hold on;
end

% ==================================================================
% for the third scalar (feedback), project out the effect of the decision interval
% ==================================================================

feedbackscalars = squeeze(nanmean(timelock(4).lock.trial(:, pupilchan, find(timelock(4).lock.time > 0) ), 3));
trialinfo(:,3)  = projectout(feedbackscalars, trialinfo(:, 2));

% ==================================================================
% lastly, take the final 100 ms of each trial
% ==================================================================

endoftrlscalar = squeeze(nanmean(timelock(4).lock.trial(:, pupilchan, ...
    end-(0.1*data.fsample):end), 3));
trialinfo(:, 4) = endoftrlscalar;

% check this all went well
assert(~any(isnan(trialinfo(:))));

end

