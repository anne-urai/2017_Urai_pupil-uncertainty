function [] = fig4c_decisionStrategies(lagGroups, whichmodulator)

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'rt'; end

% ========================================================= %
% panel B: decision strategy for lags 1-3
% ========================================================= %

clc;
subjects = 1:27;
colors = cbrewer('qual', 'Set1', 9);
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));

posRespSj = find(dat.response(:, 1) > 0);
negRespSj = find(dat.response(:, 1) < 0);

%subplot(442);
hold on;

% start with text
fz = 6;
text(0, 0.36, 'win stay', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, 0.32, 'lose switch', 'horizontalalignment', 'center', 'fontsize', fz);

text(0, -0.32, 'win switch', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, -0.36, 'lose stay', 'horizontalalignment', 'center', 'fontsize', fz);

text(0.36, .05, 'stay', 'rotation', 270, 'fontsize', fz);
text(-0.36, -.06, 'switch', 'rotation', 90, 'fontsize', fz);

plot([-1 1], [-1 1], 'color', 'k', 'linewidth', 0.5);
plot([-1 1], [1 -1], 'color', 'k', 'linewidth', 0.5);

plot(mean(dat.response(posRespSj, lagGroups), 2), mean(dat.stimulus(posRespSj, lagGroups), 2), ...
    '.', 'MarkerFaceColor', colors(2, :), 'MarkerEdgeColor', colors(2, :), 'markersize', 8);
plot(mean(dat.response(negRespSj, lagGroups), 2), mean(dat.stimulus(negRespSj, lagGroups), 2), ...
    '.', 'MarkerFaceColor', colors(5, :), 'MarkerEdgeColor', colors(5, :), 'markersize', 8);
axis tight; axis square; xlims = get(gca, 'xlim'); ylims = get(gca, 'ylim');
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));


% also show the mean
h = ploterr(mean(mean(dat.response(:, lagGroups), 2)), mean(mean(dat.stimulus(:, lagGroups), 2)), ...
    std(mean(dat.response(:, lagGroups), 2)), ...
    std(mean(dat.stimulus(:, lagGroups), 2)), ...
    'o');
set(h(1), 'markeredgecolor', 'k', 'markerfacecolor', 'w', 'markersize', 4);
set(h(2), 'color', 'k');
set(h(3), 'color', 'k');

maxlim = 0.4;
xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
set(gca, 'xtick', -maxlim:maxlim:maxlim, 'ytick', -maxlim:maxlim:maxlim);
xlabel('Response weight'); ylabel('Stimulus weight');
box on; axis square;


print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4c_decisionStrategies.pdf'));
