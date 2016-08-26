close all; figure;
global mypath;

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

plot(newx, newy); hold on; plot(plotx, mean(grandavg.acc));
grid on; waitforbuttonpress;

% find the point at which y reaches 100
[~, thisy] = min(abs(newy-1));
disp(thisy);
nrBins     = floor(newx(thisy));

% now use this range to plot
newx = linspace(-nrBins, nrBins, 1000);
newy = p(1)*newx + p(2);

subplot(441);
plot(newx, newy*100, 'k'); 
hold on; 
h = boundedline(plotx, 100* nanmean(grandavg.acc), ...
    100 * nanstd(grandavg.acc) ./ sqrt(length(subjects)), 'cmap', [0 0 0]);

axis square; box off;
ylim([45 100]);
plot([min(newx) max(newx)], [50 50], 'color', [0.5 0.5 0.5]);
xlim([min(newx)*1.1 max(newx)]);
set(gca, 'ytick', [50 75 100], 'xtick', [min(newx) min(plotx) max(plotx) max(newx)]);
set(gca, 'xticklabel', {'lowest', 'low', 'high', 'highest'});
xlabel('Pupil response');
ylabel('Accuracy (%)');

print(gcf, '-dpdf', sprintf('%s/Figures/figureS3.pdf', mypath));

