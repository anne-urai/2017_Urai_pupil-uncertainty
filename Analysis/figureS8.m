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

clearvars -except mypath;
global mypath
close all; clc;

nbins = 3;
grandavg = postPupilBehaviour('baseline_pupil', nbins, []);
subplot(4,4,1);
% line to indicate 50% repetition
y = grandavg.repetition;
plot([1 nbins], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;
thiscolor = [0 0 0]; thismarker = '.'; thismarkersize = 14;
% errorbar
h = ploterr(1:nbins, nanmean(y), [], nanstd(y) ./sqrt(27), 'k-',  'abshhxy', 0);
set(h(1), 'color', thiscolor, 'markersize', thismarkersize, 'marker', thismarker);
set(h(2), 'color', thiscolor); % line color

xticklabs       = repmat({' '}, 1, nbins);
xticklabs{1}    = 'low';
xticklabs{end}  = 'high';
if nbins == 3, xticklabs{2} = 'med'; end

set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
    'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
% axis square;
xlim([0.5 nbins+0.5]);

% determine y label and limits
set(gca, 'ylim', [0.45 0.55], 'ytick', 0.45:0.05:0.55);
ylabel('P(repeat)');

% do Bayesian ANOVA to get Bayes Factors
statdat             = table;
statdat.DV          = y(:);
sj = repmat(1:27, nbins, 1)';
statdat.subjnr      = sj(:);
ft      = repmat(transpose(1:nbins), 1, 27)';
statdat.prevPupilBins = ft(:);
writetable(statdat, sprintf('%s/Data/CSV/ANOVAdat.csv', mypath));
system('/Library/Frameworks/R.framework/Resources/bin/R < BayesFactorANOVA.R --no-save');
statres = readtable(sprintf('%s/Data/CSV/ANOVAresults.csv', mypath)); % fetch results

yval = max(get(gca, 'ylim'));
if statres.pvalue < 0.05, % only show if significant
    mysigstar(gca, [1 nbins], [yval yval], statres{cnt}.pvalue, 'k', 'down');
end
fprintf('\n ANOVA F(%d,%d) = %.3f, p = %.3f, bf10 = %.3f \n', ...
     statres.df1, statres.df2, statres.F, statres.pvalue, statres.bf10);
xlabel({'Current baseline'; 'pupil diameter'}); axis square;

%% all 7 lags for pupil x choice weight
subplot(442);
plot([1 7], [0 0], ':'); hold on;
set(gca, 'ycolor', 'k', 'xtick', 1:7, 'ytick', [-0.1:0.05:0.1], ...
    'ylim', [-.07 .01], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
ylabel('Pupil x choice weight');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
xlabel('Lags'); axis square;
plot(1, -0.06, '*k', 'markersize', 2);
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
boundedline(1:7, mean(dat.response_pupil), std(dat.response_pupil) ./ sqrt(27), 'k');

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS8.pdf', mypath));
