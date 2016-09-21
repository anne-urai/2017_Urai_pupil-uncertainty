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

close all; figure;
global mypath;

% error vs correct
subplot(441); uncertainty_byErrorCorrect('decision_pupil');
% subplot(4,7,10); plotBetasSwarm(b, colors([1 2], :));
% set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'});


% psychometric functions
subplot(4,4,2);  psychFuncs_byUncertainty('decision_pupil');
%subplot(4,7,24); plotBetasSwarm(1./b, [0.7 0.7 0.7; 0.2 0.2 0.2]);
%set(gca, 'xtick', [1 2], 'xticklabel', {'low', 'high'});
%xlabel('Pupil response'); ylabel('Sensitivity (a.u.)');
%ylim([0 1]);

% other metrics of uncertainty
subplot(4,4,3); uncertaintyAccuracy('decision_pupil');
% subplot(4,7,17); plotBetasSwarm(b(:, 2), [0 0 0]);
% set(gca, 'xtick', 1, 'xticklabel', []);

% INSET, HAS TO BE MADE SMALLER IN ILLUSTRATOR
subplot(444);
nbins = 9;
plotx = -floor(nbins/2) : floor(nbins/2);
% can try this also with all subjects
subjects = 1:27;
field           = 'decision_pupil';
grandavg.pup    = nan(length(subjects), nbins);
grandavg.acc    = nan(length(subjects), nbins);
grandavg.b      = nan(length(subjects), 2);

for sj = subjects,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % bin
    [grandavg.pup(sj, :), grandavg.acc(sj, :)] = divideintobins(data.(field), data.correct, nbins);
    
    p = polyfit(plotx, (grandavg.acc(sj, :)), 1);
    newx = linspace(-nbins*5, nbins*5, 1000);
    newy = p(1)*newx + p(2);

    % find the point at which y reaches 100
    [~, thisy] = min(abs(newy-1));
    disp(thisy);
    nrBins     = floor(newx(thisy));
    
    % now use this range to plot
    newx = linspace(-nrBins, nrBins, 1000);
    newy = p(1)*newx + p(2);
    grandavg.perf50(sj) = min(newy);
    
    % ============================================= %
    % project RT out of the pupil and vice versa
    % ============================================= %
    
    switch field
        case 'rt'
            data.(field) = projectout(zscore(data.rtNorm), zscore((data.decision_pupil)));
        case 'decision_pupil'
            data.(field) = projectout(zscore(data.decision_pupil), zscore(data.rtNorm));
    end
    
    % also do logistic regression
    grandavg.b(sj, :) = glmfit(zscore(data.(field)), data.correct, ...
        'binomial','link','logit');
end

% ============================================= %
% EXTEND TO A RANGE OF 100 TO 50 %
% ============================================= %

p = polyfit(plotx, mean(grandavg.acc), 1);
newx = linspace(-nbins*5, nbins*5, 1000);
newy = p(1)*newx + p(2);

% find the point at which y reaches 100
[~, thisy] = min(abs(newy-1));
nrBins     = floor(newx(thisy));

% now use this range to plot
newx = linspace(-nrBins, nrBins, 1000);
newy = p(1)*newx + p(2);
colors = cbrewer('seq', 'Greens', 3);
plot(newx, newy*100, '-.', 'color', colors(3, :)); 
hold on; 
% add the
h = boundedline(plotx, 100* nanmean(grandavg.acc), ...
    100 * nanstd(grandavg.acc) ./ sqrt(length(subjects)), 'cmap', [0 0 0]);

axis square; box off;
ylim([45 100]);
plot([min(newx) max(newx)], [50 50], 'color', [0.5 0.5 0.5]);
xlim([min(newx)*1.1 max(newx)]);
set(gca, 'ytick', [50 75 100], 'xtick', [min(plotx) max(plotx)]);
set(gca, 'xticklabel', {'low', 'high'});
xlabel('Pupil response');
ylabel('Accuracy (%)');

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS3.pdf', mypath));

