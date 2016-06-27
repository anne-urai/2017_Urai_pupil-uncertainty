function [] = plotBetasSwarm(beta, colors)
% using the timewindow that is indicated in the regression timecourse plot,
% show the shape of the pupil vs motionstrength pattern

if ~exist('colors', 'var'); colors = cbrewer('qual', 'Set1', 8); end
global mypath;

hold on;
for i = 1:size(beta, 2),
    barcol = (colors(i, :) * 2); % lighter
    barcol(barcol > 1) = 1;
    bar(i, squeeze(mean(beta(:, i))), 'edgecolor', 'none', ...
        'facecolor', [0.9 0.9 0.9], 'barwidth', 0.5);
    scatter(i * ones(1, size(beta, 1)), beta(:, i), ...
        3, colors(i, :), 'jitter','on', 'jitterAmount', 0.2);
end

set(gca, 'xtick', [1 2], 'xminortick', 'off');
ylabel('Beta weights (a.u.)'); 
axis tight;
xlim([0.5 size(beta,2) + 0.5]);
if size(beta, 2) == 1,
    xlim([0 2]);
end
ylims = get(gca, 'ylim');
yrange = range(ylims);
ylim([ylims(1) - yrange*0.1 ylims(2) + yrange*0.1]);

% stats
[~, pval(1), ~, stat] = ttest(beta(:, 1), 0, 'tail', 'both');
mysigstar(gca, 1, max(get(gca, 'ylim')), pval(1));

if size(beta,2) == 2,
    [~, pval(2), ~, stat] = ttest(beta(:, 2), 0, 'tail', 'both');
    [~, pval(3), ~, stat] = ttest(beta(:,1), beta(:,2));
    
    mysigstar(gca, 2, max(get(gca, 'ylim')), pval(2));
    mysigstar(gca, [1 2], max(get(gca, 'ylim'))*1.1, pval(3));
end

end
