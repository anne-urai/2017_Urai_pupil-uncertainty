function [] = fig4b_decisionStrategies(whichmodulator)
% decision strategy space

global mypath;

hold on;
plot([-1 1], [-1 1], 'color', 'k', 'linewidth', 0.5);
plot([-1 1], [1 -1], 'color', 'k', 'linewidth', 0.5);

load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator));

for sj = 1:27,
    h = ploterr(dat.response(sj, 1), dat.stimulus(sj, 1), ...
        {dat.responseCI(sj, 1, 1) dat.responseCI(sj, 1, 2)}, ...
        {dat.stimulusCI(sj, 1, 1) dat.stimulusCI(sj, 1, 2)}, '.', 'abshhxy', 0);
    set(h(1), 'color', mycolmap(sj, :), 'markerfacecolor', mycolmap(sj, :));
    set(h(2), 'color', mycolmap(sj, :), 'linewidth', 0.5);
    set(h(3), 'color', mycolmap(sj, :), 'linewidth', 0.5);
end

% also show the mean
h = ploterr(mean(mean(dat.response(:, 1), 2)), mean(mean(dat.stimulus(:, 1), 2)), ...
    std(mean(dat.response(:, 1), 2)), ...
    std(mean(dat.stimulus(:, 1), 2)), ...
    'o', 'abshhxy', 0);
set(h(1), 'markeredgecolor', 'k', 'markerfacecolor', 'w', 'markersize', 4);
set(h(2), 'color', 'k');
set(h(3), 'color', 'k');

% start with text
fz = 6;
text(0, 0.36, 'win stay', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, 0.32, 'lose switch', 'horizontalalignment', 'center', 'fontsize', fz);

text(0, -0.32, 'win switch', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, -0.36, 'lose stay', 'horizontalalignment', 'center', 'fontsize', fz);

text(0.36, .05, 'stay', 'rotation', 270, 'fontsize', fz);
text(-0.36, -.06, 'switch', 'rotation', 90, 'fontsize', fz);

% layout
maxlim = 0.4;
xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
maxlim = 0.4;
set(gca, 'xtick', -maxlim:maxlim:maxlim, 'ytick', -maxlim:maxlim:maxlim);
xlabel('Choice weight'); ylabel('Stimulus weight');
box on; axis square;

end