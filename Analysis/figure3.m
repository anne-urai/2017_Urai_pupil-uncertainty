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

close all; 
global mypath;

subplot(4,4,1); pupilTimecourse(0);
subplot(4,4,3); pupilUncertaintyTimecourse;

% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 9);
% error vs correct
subplot(445); b = uncertainty_byErrorCorrect('decision_pupil');
cla; plotBetasSwarm(b, colors([1 2], :));
set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'});

% other metrics of uncertainty
subplot(4,4,6); b = uncertaintyAccuracy('decision_pupil');
cla; plotBetasSwarm(b(:, 2), [0 0 0]);
set(gca, 'xtick', 1, 'xticklabel', []);

% psychometric functions
subplot(4,4,7); b = psychFuncs_byUncertainty('decision_pupil');
cla; plotBetasSwarm(1./b, [0.7 0.7 0.7; 0.2 0.2 0.2]);
set(gca, 'xtick', [1 2], 'xticklabel', {'low', 'high'});
xlabel('Pupil response'); ylabel('Sensitivity (a.u.)');
ylim([0 1]);

% ensure same axes proportions
ax = findobj(gcf, 'type', 'axes');
for a = 1:length(ax),
    pos = get(ax(a), 'position'); pos(4) = 0.15; set(ax(a), 'position', pos);
end
print(gcf, '-dpdf', sprintf('%s/Figures/Figure3.pdf', mypath));
% note: several of these panels appear in Figure S3 in the paper.

%% compute the correlation between RT and pupil values for each SJ
clc;
subjects = 1:27; 
for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    [thispearson(sj), thispearsonpval(sj)] = corr(data.decision_pupil, data.rtNorm, 'type', 'pearson', 'rows', 'complete');
end
% check across the group
fprintf('average r: %.3f range: %.3f to %.3f, significant in %d out of %d observers \n', ...
    mean(thispearson), min(thispearson), max(thispearson), length(find(thispearsonpval < 0.05)), 27)
