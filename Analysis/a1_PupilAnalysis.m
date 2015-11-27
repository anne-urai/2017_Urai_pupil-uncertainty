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

% add paths and fieldtrip
close all; clc;
addpath(genpath('~/code/pupilUncertainty'));
addpath(genpath('~/code/Tools'));
addpath(genpath('~/Dropbox/code/pupilUncertainty'));
addpath(genpath('~/Dropbox/code/Tools'));
addpath('~/Documents/fieldtrip/');
ft_defaults;

try
    pathname = '~/Data/UvA_pupil';
    cd(pathname);
catch % on the cluster
    pathname = '~/Data/HD1/UvA_pupil';
    cd(pathname);
end

% check if this file doesnt exist yet
if exist(sprintf('%s/P%02d_alleye2.mat', pathname, sj), 'file');
 %   disp('skipping');
 %   return
end

% subject specific folder call P01, with one session S1-S6 containing all
% the pupil files
cd(sprintf('%s/P%02d/', pathname, sj));

clear sessions;
% check which sessions to use
s = dir('S*');
s = {s(:).name};
for i = 1:length(s), sessions(i) = str2num(s{i}(2)); end

for session = unique(sessions),
    
    cd([sprintf('%s/P%02d/', pathname, sj) 'S' num2str(session)]);
    
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
        clearvars -except sj session block subjects sessions blocks pathname
        
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
        data = blink_regressout(data, blinksmp, saccsmp, 1);
        saveas(gcf,  sprintf('%s/Figures/P%02d_s%d_b%d_projectout.pdf', pathname, sj, session, block), 'pdf');

        % ==================================================================
        % compute percent signal change rather than zscore
        % median is less sensitive to outliers
        % ==================================================================
        
        pupildat    = data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:);
        medianpupil = median(pupildat);
        
        % normalize
        pupiltimecourse = (pupildat - medianpupil) ./ medianpupil * 100; % normalize
        data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:) = pupiltimecourse; % put back in
        
        % ==================================================================
        % define trials
        % ==================================================================
        
        cfg                         = [];
        cfg.trialfun                = 'trialfun_EL_UvA';
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
        disp(['Saving... ' sprintf('P%02d_s%d_b%02d_eyeclean2.mat', sj, session, block)]);
        % save these datafiles before appending
        savefast(sprintf('P%02d_s%d_b%02d_eyeclean2.mat', sj, session, block), 'data');
        cd(['S' num2str(session)]);
        
        %  end
    end
end

% ==================================================================
% now append all the eyelink files together
% ==================================================================

% check if the full dataset is not there yet
cd(sprintf('%s/P%02d/', pathname, sj));
eyelinkfiles = dir(sprintf('P%02d*_eyeclean2.mat', sj));

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
cfg.outputfile = sprintf('%s/P%02d_alleye2.mat', pathname, sj);
ft_appenddata(cfg);

end
