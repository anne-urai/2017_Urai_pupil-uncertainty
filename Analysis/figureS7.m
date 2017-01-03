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

close;
figure;

% combine the 2 models, where pupil and RT were separately added as
% interaction variables, to redo figure 5
dat1 = load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil'));
dat2 = load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'rt'));
dat = dat1.dat; 
dat.response_rt = dat2.dat.response_pupil;
dat.stimulus_rt = dat2.dat.stimulus_pupil;
dat.response_rtCI = nan(size(dat1.dat.response_pupilCI));
dat.response_pupilCI = nan(size(dat1.dat.response_pupilCI));
dat.responseCI = zeros(size(dat.responseCI));
save(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt_separate'));

subplot(4,4,9); plotBetas([dat.response_pupil(:, 1) ...
    dat.stimulus_pupil(:, 1)  dat.response_rt(:, 1) dat.stimulus_rt(:, 1)], ...
    0.5*ones(4,3));

% paired stats
[~, pval, ~, stat] = ttest(dat.response_pupil(:, 1), dat.stimulus_pupil(:, 1));
mysigstar(gca, [1 2], -0.08, pval);
[~, pval, ~, stat] = ttest(dat.response_rt(:, 1), dat.stimulus_rt(:, 1));
mysigstar(gca, [3 4], -0.08, pval);

xlim([0.5 4.5]);
set(gca, 'xtick', 1:4, 'xticklabel', ...
    {'Pupil x choice', 'Pupil x stimulus', 'RT x choice', 'RT x stimulus'}, ...
    'xticklabelrotation', -30); axis square;

subplot(4,4,10); sjCorrelation('pupil', 'response', 'pupil+rt_separate');
subplot(4,4,11); sjCorrelation('rt', 'response', 'pupil+rt_separate');

% show median split for correlation stuff
subplot(4,4, 15); medianSplit('pupil', 'response',  'pupil+rt_separate');
ylim([-0.4 0.15]); ylabel('Pupil x choice weight');
axis square;
subplot(4,4,16 ); medianSplit('rt', 'response',  'pupil+rt_separate');
ylim([-0.4 0.15]); set(gca, 'yaxislocation', 'right');
ylabel('RTx choice weight');
axis square;

print(gcf, '-dpdf',  sprintf('%s/Figures/FigureS7.pdf', mypath));
