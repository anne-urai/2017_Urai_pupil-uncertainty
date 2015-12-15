function [] = fig4c_decisionStrategies_interaction(lagGroups, whichmodulator)

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end

% ========================================================= %
% panel B: decision strategy for lags 1-3
% ========================================================= %

clc;
subjects = 1:27;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
colors = linspecer(9);

posRespSj = find(dat.response(:, 1) > 0);
negRespSj = find(dat.response(:, 1) < 0);

load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));
%subplot(442);
hold on;

plot([-1 1], [-1 1], 'color', 'k');
plot([-1 1], [1 -1], 'color', 'k');

plot(mean(dat.response_pupil(posRespSj, lagGroups), 2), mean(dat.stimulus_pupil(posRespSj, lagGroups), 2), ...
    'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(8, :), 'markersize', 3);
plot(mean(dat.response_pupil(negRespSj, lagGroups), 2), mean(dat.stimulus_pupil(negRespSj, lagGroups), 2), ...
    'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(9, :), 'markersize', 3);
axis tight; axis square; xlims = get(gca, 'xlim'); ylims = get(gca, 'ylim');

% also show the mean
h = ploterr(mean(mean(dat.response_pupil(:, lagGroups), 2)), mean(mean(dat.stimulus_pupil(:, lagGroups), 2)), ...
    std(mean(dat.response_pupil(:, lagGroups), 2)), ...
    std(mean(dat.stimulus_pupil(:, lagGroups), 2)), ...
    'o',  'hhxy', 0.0000001);
set(h(1), 'markeredgecolor', colors(2,:), 'markerfacecolor', 'w', 'markersize', 4);
set(h(2), 'color', colors(2,:));
set(h(3), 'color', colors(2,:));

maxlim = 0.12;
xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
set(gca, 'xtick', -maxlim:maxlim:maxlim, 'ytick', -maxlim:maxlim:maxlim);
xlabel('Pupil * response'); ylabel('Pupil * stimulus');
box on; axis square;

if 0,
    fz = 6;
    text(0, 0.26, 'win stay', 'horizontalalignment', 'center', 'fontsize', fz);
    text(0, 0.22, 'lose switch', 'horizontalalignment', 'center', 'fontsize', fz);
    
    text(0, -0.22, 'win switch', 'horizontalalignment', 'center', 'fontsize', fz);
    text(0, -0.26, 'lose stay', 'horizontalalignment', 'center', 'fontsize', fz);
    
    text(0.26, .05, 'stay', 'rotation', 270, 'fontsize', fz);
    text(-0.26, -.06, 'switch', 'rotation', 90, 'fontsize', fz);
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4c_decisionStrategiesInteraction.pdf'));
