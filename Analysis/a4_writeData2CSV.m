function [] = a4_writeData2CSV()
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

% from the pupil file and the motionenergy overview, make one csv file to
% quickly do behavioural analyses per subject
% the pupil data will be segmented here, and a single-trial scalar for the
% pupil response is defined based on the significant time window obtained
% in fig3b.

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
    earlyDecisionPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, :), 3));
    
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
        'decision_pupil', 'feedback_pupil', 'trialend_pupil', 'response_pupil'});

    writetable(t, sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    disp(['finished sj ' num2str(sj)]);
    alldat{find(sj==subjects)} = t;
    toc;
end

% also make 1 big file for all subjects
disp('writing to file...');
alldat2 = cat(1, alldat{:});
writetable(alldat2, sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));

end
