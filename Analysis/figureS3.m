close all; figure;
global mypath;

% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 9);
% error vs correct
subplot(441); Uncertainty_byErrorCorrect('decision_pupil');
% subplot(4,7,10); plotBetasSwarm(b, colors([1 2], :));
% set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'});


% psychometric functions
subplot(4,4,2);  PsychFuncs_byUncertainty('decision_pupil');
%subplot(4,7,24); plotBetasSwarm(1./b, [0.7 0.7 0.7; 0.2 0.2 0.2]);
%set(gca, 'xtick', [1 2], 'xticklabel', {'low', 'high'});
%xlabel('Pupil response'); ylabel('Sensitivity (a.u.)');
%ylim([0 1]);

% other metrics of uncertainty
subplot(4,4,3); UncertaintyAccuracy('decision_pupil');
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

