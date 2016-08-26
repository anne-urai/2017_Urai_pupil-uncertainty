function [] = a4_writeData2CSV()
% from the pupil file and the motionenergy overview, make one csv file to
% quickly do behavioural analyses per subject
% the pupil data will be segmented here, and a single-trial scalar for the
% pupil response is defined based on the significant time window obtained
% in fig3b.
%
% Anne Urai, 2015
global mypath;

% we already have all the info we need in the grandaverage file
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));

subjects = 1:length(pupilgrandavg.timelock);
for sj = (subjects),
    tic;
    
    clearvars -except sj subjects alldat pupilgrandavg mypath;
    
    % ==================================================================
    % get trialinfo as it is
    % ==================================================================
    
    newtrl = pupilgrandavg.timelock{sj}(4).lock.trialinfo;
    newtrl = [newtrl(:, 1:11) sj*ones(length(newtrl), 1) newtrl(:, 12:15)];
    
    % ==================================================================
    % get pupil scalars
    % ==================================================================
    
    data.fsample    = 100; % make sure we use the resampled frequency from the pupilAnalysis
    pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(4).lock.label, 'EyePupil')==1);
    
    % around the decision peak
    earlyDecisionPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, ...
        find(pupilgrandavg.timelock{sj}(3).lock.time > 0 & pupilgrandavg.timelock{sj}(3).lock.time < 1 ) ), 3));
    
    % 250 ms before feedback
    decisionPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, ...
        find(pupilgrandavg.timelock{sj}(4).lock.time < 0 & pupilgrandavg.timelock{sj}(4).lock.time > -0.250 ) ), 3));
    
    % peak of error vs correct betas across the group is at 640 ms after feedback
    % so for the scalar, we take 515 - 765ms after fb
    feedbackPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, ...
        find(pupilgrandavg.timelock{sj}(4).lock.time > 0.515 & pupilgrandavg.timelock{sj}(4).lock.time < 0.765 ) ), 3));
    
    endofTrlPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, ...
        end-(0.1*data.fsample):end), 3));
    
    % put together
    newtrl = [newtrl decisionPupil feedbackPupil endofTrlPupil earlyDecisionPupil];
    
    % ==================================================================
    % write to csv for each individual
    % ==================================================================
    
    t = array2table(newtrl, 'VariableNames', ...
        {'stim', 'coherence',  'difficulty', 'motionstrength', ...
        'resp', 'rt', 'correct', 'correctM', ...
        'trialnr', 'blocknr', 'sessionnr', 'subjnr',  ...
        'baseline_pupil', ...
        'int1motion', 'int2motion', 'rtNorm', ...
        'decision_pupil', 'feedback_pupil', 'trialend_pupil', 'earlydecision_pupil'});

    writetable(t, sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    disp(['finished sj ' num2str(sj)]);
    alldat{find(sj==subjects)} = t;
    toc;
end

if length(subjects) > 5,
    disp('writing to file...');
    alldat2 = cat(1, alldat{:});
    
    %     % write to csv for all subjects
    %     t = array2table(alldat2, 'VariableNames', ...
    %         {'stim', 'coherence',  'difficulty', 'motionstrength', ...
    %         'resp', 'rt', 'correct', 'correctM', ...
    %         'trialnr', 'blocknr', 'sessionnr', 'subjnr',  ...
    %         'baseline_pupil', 'decision_pupil', 'feedback_pupil', 'trialend_pupil', 'earlydecision_pupil'});

    writetable(alldat2, sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
end

end
