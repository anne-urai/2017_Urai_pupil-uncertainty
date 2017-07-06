function [] = pupilTimecourse_feedbackLocked()
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

global mypath;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPLIT BY CORR VS ERROR AND DIFFICULTY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except mypath plotTotalAverage; subjects = 1:27;
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens

% get all data
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));

% append all the mean timecourses per condition
for sj = unique(subjects),
    
    pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'EyePupil')==1);
    
    % baseline correct at feedback
    decisionPupil = squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, ...
        find(pupilgrandavg.timelock{sj}(4).lock.time < 0 & pupilgrandavg.timelock{sj}(4).lock.time > -0.5) ), 3));
    pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :) = bsxfun(@minus, ...
        pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :), decisionPupil);
    
    thistable = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    cors = [0 1];
    cnt  = 0;
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
            alltimelock(sj, cnt, :) = squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(trls, pupilchan, :)));
        end
    end
end

% color scheme
cols = cbrewer('div', 'RdYlBu', 10);
cols = cols([3 2 1 end-2:end], :);

% plot
hold on;
plot([0 0], [-0.5 0.5], 'color', [0.5 0.5 0.5]);
plot([-1 4], [0 0], 'color', [0.5 0.5 0.5]);

ph = boundedline(pupilgrandavg.timelock{end}(4).lock.time, squeeze(nanmean(alltimelock)), ...
    permute(squeeze(nanstd(alltimelock)) / sqrt(length(subjects)), [2 3 1]), ...
    'cmap', cols);
ylabel('Pupil response (z)');
axis tight;
legnames = {'error weak', 'error medium', 'error strong', 'correct weak', 'correct medium', 'correct strong'};
lh = legend(ph, legnames); % make t
lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .25;
set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 8);
set(gca, 'ytick', [-.5 0 0.5], 'ylim', [-0.6 0.5], 'xlim', [-1.1 4]);
xlabel('Time from feedback (s)');
set(gca, 'xtick', -1:1:4);

end

