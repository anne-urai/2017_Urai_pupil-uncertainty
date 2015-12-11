function [] = a4_writeData2CSV()
% from the pupil file and the motionenergy overview, make one csv file to
% quickly do behavioural analyses per subject
% the pupil data will be segmented here, and a single-trial scalar for the
% pupil response is defined based on the significant time window obtained
% in fig3b.
%
% Anne Urai, 2015

close all; 
addpath('~/Dropbox/code/pupilUncertainty/Analysis');
dbstop if error;
addpath('~/Documents/fieldtrip');
ft_defaults;

% we already have all the info we need in the grandaverage file
load('~/Data/pupilUncertainty/GrandAverage/pupilgrandaverage.mat');

subjects = 1:length(pupilgrandavg.timelock);
for sj = fliplr(subjects),
    tic;
    
    clearvars -except sj subjects alldat pupilgrandavg;
    
    % ==================================================================
    % get trialinfo as it is
    % ==================================================================
    
    newtrl = pupilgrandavg.timelock{sj}(4).lock.trialinfo;
    newtrl = [newtrl(:, 1:11) sj*ones(length(newtrl), 1) newtrl(:, 12)];
    
    % ==================================================================
    % get pupil scalars
    % ==================================================================
    
    data.fsample          = 100; % make sure we use the resampled frequency from the pupilAnalysis
   
    % baseline, decision, feedback, end of trial
    baselinePupil = pupilgrandavg.timelock{sj}(1).lock.bl';
    
    pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(4).lock.label, 'EyePupil')==1);

    % 250 ms before feedback
    decisionPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, ...
        find(pupilgrandavg.timelock{sj}(4).lock.time < 0 & pupilgrandavg.timelock{sj}(4).lock.time > -0.250 ) ), 3));
    
    % 250 ms after feedback
    feedbackPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, ...
        find(pupilgrandavg.timelock{sj}(4).lock.time > 0 & pupilgrandavg.timelock{sj}(4).lock.time < 0.250 ) ), 3));
    
    endofTrlPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, ...
    end-(0.1*data.fsample):end), 3));

    % put together
    newtrl = [newtrl decisionPupil feedbackPupil endofTrlPupil];
    
    % ==================================================================
    % write to csv for each individual
    % ==================================================================
    
    t = array2table(newtrl, 'VariableNames', ...
        {'stim', 'coherence',  'difficulty', 'motionstrength', ...
        'resp', 'rt', 'correct', 'correctM', ...
        'trialnr', 'blocknr', 'sessionnr', 'subjnr',  ...
        'baseline_pupil', 'decision_pupil', 'feedback_pupil', 'trialend_pupil'});
    
    writetable(t, sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    disp(['finished sj ' num2str(sj)]);
    alldat{find(sj==subjects)} = newtrl;
    
end

if length(subjects) > 5,
    disp('writing to file...');
    alldat2 = cat(1, alldat{:});
    
    % write to csv for all subjects
    t = array2table(alldat2, 'VariableNames', ...
        {'stim', 'coherence',  'difficulty', 'motionstrength', ...
        'resp', 'rt', 'correct', 'correctM', ...
        'trialnr', 'blocknr', 'sessionnr', 'subjnr',  ...
        'baseline_pupil', 'decision_pupil', 'feedback_pupil', 'trialend_pupil'});
    
    writetable(t, sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_allsj.csv'));
end

end
