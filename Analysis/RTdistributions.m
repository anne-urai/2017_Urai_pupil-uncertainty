function RTdistributions

global mypath;
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));

colors = cbrewer('qual', 'Set1', 8);
cols = colors([1 2], :);

histogram(data.rt(data.correct == 1), 100, ...
    'displaystyle', 'stairs', 'edgecolor', cols(2, :), 'linewidth', 1, ...
    'linestyle', '-');
hold on;
histogram(data.rt(data.correct == 0), 100, ...
    'displaystyle', 'stairs', 'edgecolor', cols(1, :), 'linewidth', 1, ...
    'linestyle', '-');

xlim([-.1 2.5]); xlabel('Reaction time (s)'); ylabel('Count (x1000)'); box off;
set(gca, 'ytick', [0 2000 4000], 'yticklabel', [0 2 4]);
set(gca, 'xminortick', 'on');
axis square;

% add text + arrow instead of legend?
legend({'correct', 'error'}, 'location', 'northeast');
legend boxoff;

end