% ========================================================= %
% panel B: decision strategy for lags 1-3
% ========================================================= %

clear all; clc; close all;
whichmodulator = 'pupil';
subjects = 1:27;
lagGroups = 1:3;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

subplot(331); hold on;

plot([-1 1], [-1 1], 'color', 'k');
plot([-1 1], [1 -1], 'color', 'k');

plot(mean(dat.response(:, lagGroups), 2), mean(dat.stimulus(:, lagGroups), 2), ...
    '.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', [0.8 0.8 0.8], 'markersize', 10);
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

text(0, 0.26, 'win stay, lose switch', 'horizontalalignment', 'center');
text(0, -0.26, 'win switch, lose stay', 'horizontalalignment', 'center');
text(0.26, .04, 'stay', 'rotation', 270);
text(-0.26, -.06, 'switch', 'rotation', 90);

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/figS_decisionStrategies.pdf'));
