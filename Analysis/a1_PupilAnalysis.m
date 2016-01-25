function a1_PupilAnalysis(sj)
% pupil preprocessing and analysis pipeline. For the participant specified,
% all edf/asc EyeLink files will be read in, blink and saccade responses
% will be removed, and data will be epoched and put into a FieldTrip-like
% structure with a trialinfo matrix that will later be used to create csv
% files.
%
% Anne Urai, 2015

% if we're running this on torque, make sure the input arg is a number
if ischar(sj), sj = str2double(sj); end

% subject specific folder call P01, with one session S1-S6 containing all
% the pupil files
cd(sprintf('%s/Data/P%02d/', mypath, sj));

% check which sessions to use
s = dir('S*');
s = {s(:).name};
for i = 1:length(s), sessions(i) = str2num(s{i}(2)); end

for session = unique(sessions),
    
    cd([sprintf('%s/Data/P%02d/', pathname, sj) 'S' num2str(session)]);
    delete *.png
    
    % ==================================================================
    % LOAD IN SUBJECT SPECIFICS AND READ DATA
    % ==================================================================
    
    blocks = 1:10;
    
    % some subjects didnt do all blocks, manually correct
    if sj == 2 && session == 1,
        blocks = 4:9;
    elseif sj == 17 && session == 1,
        blocks = 1:5;
    elseif sj == 15 && session == 3,
        blocks = 1:5;
    elseif sj == 15 && session == 6,
        blocks = 1:5;
    end
    
    for block = unique(blocks),
        clearvars -except sj session block subjects sessions blocks pathname regressout
        
        disp(['Analysing subject ' num2str(sj) ', session ' num2str(session) ', block ' num2str(block)]);
        
        edffile   = dir(sprintf('P%d_s%d_b%d_*.edf', sj, session, block));
        ascfile   = dir(sprintf('P%d_s%d_b%d_*.asc', sj, session, block));
        
        % specify the filename
        if ~exist(ascfile.name, 'file'),
            
            % IF NECESSARY, CONVERT TO ASC
            if exist('~/code/Tools/eye/edf2asc-linux', 'file'),
                system(sprintf('%s %s -input', '~/code/Tools/eye/edf2asc-linux', edffile.name));
            else
                system(sprintf('%s %s -input', '~/Dropbox/code/Tools/eye/edf2asc-mac', edffile.name));
            end
            ascfile   = dir(sprintf('P%d_s%d_b%d_*.asc', sj, session, block));
        end
        
        % ==================================================================
        % making a FieldTrip structure out of EyeLink data
        % ==================================================================
        
        clear blinksmp saccsmp
        load(sprintf('P%02d_s%d_b%02d_eye.mat', sj, session, block));
        
        if ~exist('blinksmp', 'var') || ~exist('saccsmp', 'var'),
            % read in the asc EyeLink file
            asc = read_eyelink_ascNK_AU(ascfile.name);
            
            % create events and data structure, parse asc
            [data, event, blinksmp, saccsmp] = asc2dat(asc);
            
            % save
            savefast(sprintf('P%02d_s%d_b%02d_eye.mat', sj, session, block), 'data', 'asc', 'event', 'blinksmp', 'saccsmp');
        end
        
        % ==================================================================
        % blink interpolation
        % ==================================================================
        
        newpupil = blink_interpolate(data, blinksmp, 1);
        data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:) = newpupil;
        
        suplabel(sprintf('P%02d_s%d_b%d_preproc.pdf', sj, session, block), 't');
        saveas(gcf,  sprintf('%s/Figures/P%02d_s%d_b%d_preproc.pdf', pathname, sj, session, block), 'pdf');
        
        % ==================================================================
        % regress out pupil response to blinks and saccades
        % ==================================================================
        
        % for this, use only EL-defined blinksamples
        % dont add back slow drift for now
        addBackSlowDrift = 0;
        
        data = blink_regressout(data, blinksmp, saccsmp, 1, addBackSlowDrift);
        saveas(gcf,  sprintf('%s/Figures/P%02d_s%d_b%d_projectout.pdf', pathname, sj, session, block), 'pdf');
        
        % ==================================================================
        % zscore since we work with the bandpassed signal
        % ==================================================================
        
        data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:) = ...
            zscore(data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:));
        
        % ==================================================================
        % make channels with blinks and saccades
        % ==================================================================
        
        data.label{4} = 'Blinks';
        data.trial{1}(4, :) = zeros(1, length(data.time{1}));
        for b = 1:length(blinksmp),
            data.trial{1}(4, blinksmp(b,1):blinksmp(b,2)) = 1;
        end
        
        data.label{5} = 'Saccades';
        data.trial{1}(5, :) = zeros(1, length(data.time{1}));
        for s= 1:length(saccsmp),
            data.trial{1}(5, saccsmp(s,1):saccsmp(s,2)) = 1;
        end
        
        % ==================================================================
        % define trials
        % ==================================================================
        
        cfg                         = [];
        cfg.trialfun                = 'trialfun_pupil';
        cfg.trialdef.pre            = 0;
        cfg.trialdef.post           = 4;
        cfg.event                   = event;
        cfg.dataset                 = ascfile.name;
        cfg.fsample                 = asc.fsample;
        cfg.sj                      = sj;
        cfg.session                 = session;
        [cfg]                       = ft_definetrial(cfg);
        
        data                        = ft_redefinetrial(cfg, data); %make trials
        data.trialinfo              = cfg.trl(:,4:end);
        
        % in sj 3 and 5, recode the block nrs
        if (sj == 5 && session == 1) || (sj == 3 && session == 1),
            data.trialinfo(:,13) = block;
        elseif sj == 15 && session == 6,
            data.trialinfo(:,13) = block;
        end
        
        % ==================================================================
        % downsample before saving
        % ==================================================================
        
        cfg             = [];
        cfg.resamplefs  = 100;
        cfg.fsample     = data.fsample;
        
        % see Niels' message on the FT mailing list
        samplerows = find(data.trialinfo(1,:)>100); %indices of the rows with sample values (and not event codes)
        data.trialinfo(:,samplerows) = round(data.trialinfo(:,samplerows) * (cfg.resamplefs/cfg.fsample));
        
        % use fieldtrip to resample
        data = ft_resampledata(cfg, data);
        
        cd ..
        disp(['Saving... ' sprintf('P%02d_s%d_b%02d_eyeclean.mat', sj, session, block)]);
        % save these datafiles before appending
        savefast(sprintf('P%02d_s%d_b%02d_eyeclean.mat', sj, session, block), 'data');
        cd(['S' num2str(session)]);
        
    end
end

% ==================================================================
% now append all the eyelink files together
% ==================================================================

% check if the full dataset is not there yet
cd(sprintf('%s/Data/P%02d/', pathname, sj));
eyelinkfiles = dir(sprintf('P%02d*_eyeclean.mat', sj));

% make sure these are in the right order!
% otherwise, indexing of trials will go awry
for f = 1:length(eyelinkfiles),
    scandat         = sscanf(eyelinkfiles(f).name, 'P%*d_s%d_b%d*.mat');
    snum(f,:)       = scandat';
end
[sorted, sortidx]   = sort(snum(:,1)); % sort by session
sorted(:,2)         = snum(sortidx, 2); % sort by block
eyelinkfiles        = eyelinkfiles(sortidx);

cfg = [];
cfg.inputfile = {eyelinkfiles.name};
cfg.outputfile = sprintf('%s/Data/P%02d_alleye.mat', pathname, sj);
ft_appenddata(cfg);

end

% ==================================================================
% function for epoching
% ==================================================================

function [trl, event] = trialfun_pupil(cfg)
% header and events are already in the asc structures
% returns a trl that should be identical to the structure obtained from MEG
% data

event   = cfg.event;
value   = {event(find(~cellfun(@isempty,strfind({event.value},'MSG')))).value};
sample  = [event(find(~cellfun(@isempty,strfind({event.value},'MSG')))).sample];

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * cfg.fsample);
posttrig =  round(cfg.trialdef.post * cfg.fsample);

trl = [];
session = cfg.session;
for j = 1:length(value), % loop through the trials and create the trial matrix on each trl
    
    try
        % check that this is really a fixation trigger
        if (~isempty(strfind(value{j}, 'blinkbreak_end')) && ~isempty(strfind(value{j+2}, 'ref'))) || ...
                (~isempty(strfind(value{j+1}, 'fix')) && ~isempty(strfind(value{j+2}, 'ref'))),
            
            trlbegin = sample(j) + pretrig;
            %trlend   = sample(j) + posttrig;
            offset   = pretrig;
            fixoffset = sample(j);
            
            % then find the trigger after that, corresponding to the reference
            if ~isempty(strfind(value{j+2}, 'ref')),
                refoffset = sample(j+2);
            else
                error('no refoffset sample found');
            end
            
            % find the trial nr and block nr, scan message
            scandat =  sscanf(value{j+2}, 'MSG %*f block%d_trial%d_ref70_*');
            blockcnt = scandat(1); trlcnt = scandat(2);
            
            % stimulus start
            if ~isempty(strfind(value{j+3}, 'stim_coh')),
                stimoffset = sample(j+3);
            else
                error('no stimoffset sample found');
            end
            
            % decode stimulus type
            coh =  sscanf(value{j+3}, 'MSG %*f block%*d_trial%*d_stim_coh%f_dir%*d_diff%*d');
            if coh < .7,
                stimtype = -1; % weaker
            elseif coh > .7,
                stimtype = 1; % stronger
            elseif coh == .7,
                stimtype = 0; % no difference
            end
            
            % decode stim strength
            stimstrength = single(abs(.7-coh));
            if stimstrength >= 0.019 && stimstrength <= 0.03,
                stimstrength = 0.025;
            end
            
            % difficulty level
            diff =  sscanf(value{j+3}, 'MSG %*f block%*d_trial%*d_stim_coh%*f_dir%*d_diff%d');
            if isempty(diff), diff = NaN; end
            
            % response
            if ~isempty(strfind(value{j+4}, 'resp')),
                respoffset = sample(j+4);
            else
                error('no respoffset sample found');
            end
            
            resp = sscanf(value{j+4}, 'MSG %*f block%*d_trial%*d_resp_key%f_correct%f');
            if numel(resp) == 2,
                resptype = resp(1); respcorrect = resp(2);
            elseif numel(resp) == 0,
                % some erroneous asc files without key trigger
                respcorrect = sscanf(value{j+4}, 'MSG %*f block%*d_trial%*d_resp_key%*c_correct%d');
                switch respcorrect,
                    case 1
                        resptype = stimtype;
                    case 0
                        resptype = -stimtype;
                end
            end
            
            % feedback
            if ~isempty(strfind(value{j+5}, 'feedback')),
                feedbackoffset = sample(j+5);
                
                % check feedback type
                if ~isempty(strfind(value{j+5}, 'feedback_correct1')), % correct
                    feedbacktype = 1;
                elseif ~isempty(strfind(value{j+5}, 'feedback_correct0')), % error
                    feedbacktype = 0;
                elseif ~isempty(strfind(value{j+5}, 'feedback_correctNaN')), % no response given
                    continue
                    warning('no response trial removed');
                end
            else
                warning('no feedbackoffset sample found');
                
                % load in the behavioural file that corresponds
                behavfile = dir(sprintf('~/Data/pupilUncertainty/P%02d/Behav/P%02d_s%d_*.mat', cfg.sj, cfg.sj, cfg.session));
                
                % if there are several files, continue until the right one is
                % found....
                filefound = false; cnt = 1;
                while ~filefound,
                    load(sprintf('~/Data/pupilUncertainty/P%02d/Behav/%s', cfg.sj, behavfile(cnt).name));
                    if  ~all(isnan(results.correct(blockcnt, :))),
                        filefound = true;
                    end
                    cnt = cnt + 1;
                end
                
                % get the feedbackoffset (= pupilrebound) between resp and tone
                feedbackoffset   = setup.pupilreboundtime(blockcnt, trlcnt);
                feedbackoffset   = round(respoffset + feedbackoffset*1000);
                feedbacktype     = results.correct(blockcnt, trlcnt);
            end
            
            assert(isequaln(respcorrect,feedbacktype), 'respcorrect and feedbacktype do not match');
            
            % fieldtrip allows variable trial length
            trlend = feedbackoffset + posttrig;
            
            % append all to mimic the MEG's trialinfo
            newtrl   = [trlbegin trlend offset ...
                fixoffset refoffset stimtype stimstrength diff ...
                stimoffset resptype respcorrect respoffset feedbacktype feedbackoffset trlcnt blockcnt session];
            trl      = [trl; newtrl];
        end
    catch me
        warning('aborting');
    end
end

% check that the coherence levels are properly coded
cohs = unique(trl(:,7));
if numel(cohs) ~= 5,
    error('difficulty levels not properly coded');
end

% if no difficulty was indicated in the triggers, insert based on the
% stimulus strength
for d = 1:length(cohs),
    trl(find(trl(:,7)==cohs(d)),8) = d;
end

end



