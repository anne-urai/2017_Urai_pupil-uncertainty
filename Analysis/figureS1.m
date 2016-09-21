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

close; figure;
global mypath;

edges   = 0:0.04:1.5;
distFun = @(x) histcounts(x, edges, 'normalization', 'count') ./ numel(x);
% get mean of each bin instead of edges
plotedges = mean([edges', circshift(edges, [0 1])'], 2);
plotedges = plotedges(2:end);
data    = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
data.rt(data.rt < 0) = NaN;
ylims = [-0.005 .11];

% 1. RT distribution for each subject
subplot(441); 
cmap = [0 0 0];
clear RTs; clear medians;
for sj = unique(data.subjnr)',
    RTs(sj, :) = distFun(data.rt(data.subjnr == sj));
    medians(sj) = nanmedian(data.rt(data.subjnr == sj));
end
boundedline(plotedges, squeeze(nanmedian(RTs)), permute(squeeze(iqr(RTs)) ./ sqrt(27), [2 3 1]), 'cmap', cmap);
plot([median(medians) median(medians)], [0 ylims(2)], ':', 'color', 'k');
xlabel('Response time (s)'); ylabel('Fraction of trials');
xlim([-0.1 1.5]); ylim(ylims); box off;

% 2. RT distributions for each level of evidence
subplot(442); hold on;
clear RTs; clear medians;
for sj = unique(data.subjnr)',
    for difficulty = unique(data.difficulty)',
        RTs(sj, difficulty, :)  = distFun(data.rt(data.subjnr == sj & data.difficulty == difficulty));
        medians(sj, difficulty) = nanmedian(data.rt(data.subjnr == sj & data.difficulty == difficulty));
    end
end
cmap = cbrewer('seq', 'Greens', 7); cmap = cmap(2:end, :);
h = boundedline(plotedges, squeeze(nanmedian(RTs)), permute(squeeze(iqr(RTs)) ./ sqrt(27), [2 3 1]), 'cmap', cmap);
for d = 1:5, plot([median(medians(:, d)) median(medians(:, d))], [0 ylims(2)], ':', 'color', cmap(d, :)); end
xlabel('Response time (s)'); ylabel('Fraction of trials');
xlim([-0.1 1.5]); ylim(ylims); box off;
legend(h, {'weak', '', '', '', 'strong'});

% 2. RT distributions for correct and error
subplot(443); hold on;
clear RTs; clear medians;
for sj = unique(data.subjnr)',
    cors = [0 1];
    for c = 1:2,
        RTs(sj, c, :)  = distFun(data.rt(data.subjnr == sj & data.correct == cors(c)));
        medians(sj, c) = nanmedian(data.rt(data.subjnr == sj & data.correct == cors(c)));
    end
end
cmap = cbrewer('qual', 'Set1', 2);
h = boundedline(plotedges, squeeze(median(RTs)), permute(squeeze(iqr(RTs)) ./ sqrt(27), [2 3 1]), 'cmap', cmap);
for d = 1:2, plot([median(medians(:, d)) median(medians(:, d))], [0 ylims(2)], ':', 'color', cmap(d, :)); end
xlabel('Response time (s)'); ylabel('Fraction of trials');
xlim([-0.1 1.5]); ylim(ylims); box off;
legend(h, {'error', 'corrcect'});

% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 9);
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));

% error vs correct
%subplot(445); [b, bint] = Uncertainty_byErrorCorrect('rt');
%subplot(4,8,11); plotBetasSwarm(b, colors([1 2], :));
set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'});

% other metrics of uncertainty
subplot(4,4,5); b = uncertaintyAccuracy('rt');
hold on; plot([min(get(gca, 'xtick')) max(get(gca, 'xtick'))], [50 50], ':k');
subplot(4,8,11); plotBetasSwarm(b(:, 2), [0 0 0]);
set(gca, 'xtick', 1, 'xticklabel', []);
ylim([-0.6 0]);

% psychometric functions
subplot(4,4,7); b = psychFuncs_byUncertainty('rt');
subplot(4,8,15); plotBetasSwarm(1./ b, [0.5 0.5 0.5; 0.2 0.2 0.2]);
set(gca, 'xtick', [1 2], 'xticklabel', {'fast', 'slow'});
xlabel('Reaction time'); ylabel('Sensitivity (a.u.)');
ylim([0 2]);

% ensure same axes proportions
ax = findobj(gcf, 'type', 'axes');
for a = 1:length(ax),
    pos = get(ax(a), 'position'); pos(4) = 0.15; set(ax(a), 'position', pos);
end

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS1.pdf', mypath));

%% check how many trials are faster than 200 ms
for sj = unique(data.subjnr)',
        fastResps(sj) = sum(data.rt(data.subjnr == sj) < 0.2) ...
            ./ sum(~isnan(data.rt(data.subjnr == sj)));
end

%% for rebuttal letter: cumulative distribution
% this is not in the paper
distFun = @(x) histcounts(x, edges, 'normalization', 'cdf');
subplot(4,4,16); 
cmap = [0 0 0];
clear RTs; clear medians;
for sj = unique(data.subjnr)',
    RTs(sj, :) = distFun(data.rt(data.subjnr == sj));
end
boundedline(plotedges, squeeze(nanmedian(RTs)), permute(squeeze(iqr(RTs)) ./ sqrt(27), [2 3 1]), 'cmap', cmap);
xlabel('Response time (s)'); ylabel('Cumulative probability');
xlim([-0.1 1.5]); ylim([-0.05 1]); box off;
