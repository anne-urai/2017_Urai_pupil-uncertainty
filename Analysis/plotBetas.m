function [] = plotBetas(beta, colors)
% using the timewindow that is indicated in the regression timecourse plot,
% show the shape of the pupil vs motionstrength pattern

if ~exist('colors', 'var'); colors = cbrewer('qual', 'Set1', 8); end
global mypath;

% barplot with individual datapoints
hold on;
for i = 1:size(beta, 2),
    bar(i, squeeze(nanmean(beta(:, i))), ...
        'edgecolor', 'none', 'facecolor', [0.8 0.8 0.8], 'barwidth', 0.4);
end

% add paired datapoints
% hold on;
% plot(beta', '.-', 'color', [0.6 0.6 0.6], 'linewidth', 0.2);
% axis tight; xlim([0.5 size(beta,2) + 0.5]);
% ylims = get(gca, 'ylim');
% yrange = range(ylims);
% ylim([ylims(1) - yrange*0.1 ylims(2) + yrange*0.1]);

% add error bars for SEM
h = ploterr(1:size(beta, 2), squeeze(nanmean(beta)), [], ...
    squeeze(nanstd(beta)) ./ sqrt(size(beta,1)), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none');

set(gca, 'xtick', [1 2], 'xminortick', 'off');
ylabel('Beta weights (a.u.)');

% stats
% slopes
for i = 1:size(beta, 2),
    [~, pval, ~, stat] = ttest(beta(:, i), 0, 'tail', 'both');
    mysigstar(gca, i, max(get(gca, 'ylim')), pval);
end

if size(beta,2) == 2,
    [~, pval(2), ~, stat] = ttest(beta(:, 2), 0, 'tail', 'both');
    [~, pval(3), ~, stat] = ttest(beta(:,1), beta(:,2));
    
    mysigstar(gca, 2, max(get(gca, 'ylim')), pval(2));
    mysigstar(gca, [1 2], max(get(gca, 'ylim'))*1.1, pval(3));
end

end
