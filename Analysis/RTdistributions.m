function RTdistributions

global mypath;

% for every participant, compute the distribution as a fraction of choices
for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    edges = 0:0.04:2.5;
    grandavg.corr(sj, :) = histcounts(data.rt(data.correct == 1), edges) ./height(data) * 100;
    grandavg.incorr(sj, :) = histcounts(data.rt(data.correct == 0), edges) ./height(data) * 100;
end

colors = cbrewer('qual', 'Set1', 8);
colors = colors([1 2], :); % red and blue

b = boundedline(edges(1:end-1), squeeze(mean(grandavg.incorr)), squeeze(nanstd(grandavg.incorr)) ./ sqrt(27), ...
edges(1:end-1), squeeze(mean(grandavg.corr)), squeeze(nanstd(grandavg.corr)) ./ sqrt(27), 'cmap', colors );

 xlabel('Reaction time (s)'); ylabel('% of trials'); box off;
set(gca, 'xminortick', 'on');
axis square; axis tight; ylim([0 max(get(gca, 'ylim'))])
xlim([-.1 2.5]);

% add text + arrow instead of legend?
l = legend(b, {'error', 'correct'}, 'location', 'northeast');
lpos = get(l, 'position');
lpos(1) = lpos(1) + 0.05;
set(l, 'position', lpos);
legend boxoff;

axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
  axes(a).FontSize = 7;
end

end