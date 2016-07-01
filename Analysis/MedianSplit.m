function [] = MedianSplit(whichmodulator)
% plot correlation between subjects

global mypath;
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
mycolmap = cbrewer('div', 'PuOr', 3);

% split between repeaters and alternators
repeaters = find(dat.response(:, 1) > 0);
switchers = find(dat.response(:, 1) < 0);

% plot a bargraph of the modulation weights for modulation and choice
mat = dat.(['response_' whichmodulator])(:, 1);
yval = mean(mat)-std(mat);

% significance
[~, pvals(1)] = ttest(mat(switchers));
[~, pvals(2)] = ttest(mat(repeaters));
[~, pvals(3)] = ttest2(mat(switchers), mat(repeaters));

% first, repeaters
hold on;
bar(1, mean(mat(switchers)), 'barwidth', 0.5', 'facecolor', mycolmap(1,:), 'edgecolor', 'none');
h = ploterr(1, mean(mat(switchers)), [], ...
    std(mat(switchers)) ./ sqrt(length(switchers)), 'k', 'abshhxy', 0);
set(h(1), 'marker', 'none');
mysigstar(gca, 1, yval, pvals(1));

% switchers
hold on;
bar(2, mean(mat(repeaters)), 'barwidth', 0.5', 'facecolor', mycolmap(end,:), 'edgecolor', 'none');
h = ploterr(2, mean(mat(repeaters)), [], ...
    std(mat(repeaters)) ./ sqrt(length(repeaters)), 'k', 'abshhxy', 0);
set(h(1), 'marker', 'none');
mysigstar(gca, 2, yval, pvals(2));

yval = yval*1.1;
mysigstar(gca, [1 2], [yval yval], pvals(3));

xlim([0.5 2.5]); set(gca, 'xtick', 1:2, 'xticklabel', []);
switch whichmodulator
    case 'pupil'
        ylabel('Pupil modulation weight');
    case 'rt'
        ylabel('RT modulation weight');
end
ylim([-0.12 0.02]); set(gca, 'ytick', [-0.1:0.1:0]);
axis square;

end
