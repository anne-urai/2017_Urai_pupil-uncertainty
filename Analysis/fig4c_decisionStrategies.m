function [] = fig4c_decisionStrategies(lagGroups, whichmodulator)
if ~exist('lagGroups', 'var'), lagGroups = 1:3; end

% ========================================================= %
% panel B: decision strategy for lags 1-3
% ========================================================= %

clc;
subjects = 1:27;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));
colors = linspecer(9);

posRespSj = find(mean(dat.response(:, lagGroups), 2) > 0);
negRespSj = find(mean(dat.response(:, lagGroups), 2) < 0);

subplot(443); hold on;

plot([-1 1], [-1 1], 'color', 'k');
plot([-1 1], [1 -1], 'color', 'k');

plot(mean(dat.response(posRespSj, lagGroups), 2), mean(dat.stimulus(posRespSj, lagGroups), 2), ...
    '.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', colors(8, :), 'markersize', 10);
plot(mean(dat.response(negRespSj, lagGroups), 2), mean(dat.stimulus(negRespSj, lagGroups), 2), ...
    '.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', colors(9, :), 'markersize', 10);
axis tight; axis square; xlims = get(gca, 'xlim'); ylims = get(gca, 'ylim');

% also show the mean
h = ploterr(mean(mean(dat.response(:, lagGroups), 2)), mean(mean(dat.stimulus(:, lagGroups), 2)), ...
    std(mean(dat.response(:, lagGroups), 2)), ...
    std(mean(dat.stimulus(:, lagGroups), 2)), ...
    'ko',  'hhxy', 0.0000001);
set(h(1), 'markeredgecolor', 'k', 'markerfacecolor', 'w', 'markersize', 5);

maxlim = 0.3;
xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
set(gca, 'xtick', -maxlim:maxlim:maxlim, 'ytick', -maxlim:maxlim:maxlim);
xlabel('Response weight'); ylabel('Stimulus weight');
box on; axis square;

fz = 6;
text(0, 0.26, 'win stay', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, 0.22, 'lose switch', 'horizontalalignment', 'center', 'fontsize', fz);

text(0, -0.22, 'win switch', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, -0.26, 'lose stay', 'horizontalalignment', 'center', 'fontsize', fz);

text(0.26, .05, 'stay', 'rotation', 270, 'fontsize', fz);
text(-0.26, -.06, 'switch', 'rotation', 90, 'fontsize', fz);

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4c_decisionStrategies.pdf'));
