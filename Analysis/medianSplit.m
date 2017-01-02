function [] = medianSplit(whichmodulator, field)
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
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));

% split between repeaters and alternators
mat = dat.([field '_' whichmodulator])(:, 1);

for i = 1:2,

    switch i
        case 1
            subjects{i} = find(dat.response(:, 1) < 0);
        case 2
            subjects{i} = find(dat.response(:, 1) > 0);
    end
    colors = mycolmap(subjects{i},:);
    plotBetasSwarmUnpaired(i, mat(subjects{i}, :), colors);
end

% do stats between them
xlim([0.5 2.5]);
set(gca, 'xtick', 1:2, 'xticklabel', {'alternators', 'repeaters'}, ...
   'xticklabelrotation', -30);

ylims = get(gca, 'ylim');
yrange = range(ylims);
ylim([ylims(1) - yrange*0.2 ylims(2) + yrange*0.1]);

% unpaired permutation test
pval = permtest2(mat(subjects{1}), mat(subjects{2}));
mysigstar(gca, [1 2], max(get(gca, 'ylim')), pval);

end

function [] = plotBetasSwarmUnpaired(i, beta, colors)
% using the timewindow that is indicated in the regression timecourse plot,
% show the shape of the pupil vs motionstrength pattern

if ~exist('colors', 'var'); colors = cbrewer('qual', 'Set1', 8);
colors = [0 0 0; colors]; end % start with black
hold on; % paired

bar(i, squeeze(mean(beta)), 'edgecolor', 'none', ...
    'facecolor', [0.8 0.8 0.8], 'barwidth', 0.5);

% scatter all the points
scatter(i * ones(1, size(beta, 1)), beta, ...
    10, colors, 'o', 'linewidth', 0.5, 'jitter', 'on', 'jitteramount', 0);
set(gca, 'xtick', [1 2], 'xminortick', 'off');
ylabel('Beta weights (a.u.)');

xlim([0.5 size(beta,2) + 0.5]);
axis tight;
box off;

% single sample permutation test
pval = permtest(beta);
mysigstar(gca, i, min(get(gca, 'ylim')), pval);

end
