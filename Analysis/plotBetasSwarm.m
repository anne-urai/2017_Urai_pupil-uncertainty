function [] = plotBetasSwarm(beta, colors)
% using the timewindow that is indicated in the regression timecourse plot,
% show the shape of the pupil vs motionstrength pattern

if ~exist('colors', 'var'); colors = cbrewer('qual', 'Set1', 8); 
colors = [0 0 0; colors]; end % start with black
hold on; % paired

for i = 1:size(beta, 2),
    bar(i, squeeze(mean(beta(:, i))), 'edgecolor', 'none', ...
        'facecolor', [0.8 0.8 0.8], 'barwidth', 0.5);
end

% paired lines
if size(beta, 2) == 2,
    plot(beta', '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.4);
end

% scatter all the points
for i = 1:size(beta, 2),
    scatter(i * ones(1, size(beta, 1)), beta(:, i), ...
        10, colors(i, :), 'o', 'linewidth', 0.5, 'jitter', 'on', 'jitteramount', 0);
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

box off;
end
