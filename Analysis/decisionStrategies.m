function [] = decisionStrategies(whichmodulator, errorbars, showexample)
% decision strategy space

if ~exist('errorbars', 'var'); errorbars = 1; end
if ~exist('showexample', 'var'); showexample = 1; end

global mypath;

hold on;
plot([-0.8 0.8], [-0.8 0.8], 'color', 'k', 'linewidth', 0.5);
plot([-0.8 0.8], [0.8 -0.8], 'color', 'k', 'linewidth', 0.5);

load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator));

if errorbars,
    for sj = 1:27,
        h = ploterr(dat.response(sj, 1), dat.stimulus(sj, 1), ...
            {dat.responseCI(sj, 1, 1) dat.responseCI(sj, 1, 2)}, ...
            {dat.stimulusCI(sj, 1, 1) dat.stimulusCI(sj, 1, 2)}, '.', 'abshhxy', 0);
        set(h(1), 'color', mycolmap(sj, :), 'markerfacecolor', mycolmap(sj, :));
        set(h(2), 'color', mycolmap(sj, :), 'linewidth', 0.5);
        set(h(3), 'color', mycolmap(sj, :), 'linewidth', 0.5);
    end
else
    for sj = 1:27,
        h = ploterr(dat.response(sj, 1), dat.stimulus(sj, 1), ...
            [], ...
            [], '.', 'abshhxy', 0);
        set(h(1), 'color', mycolmap(sj, :), 'markerfacecolor', mycolmap(sj, :), 'markersize', 10);
    end
end

fz = 6;
text(0, 0.36, 'win stay', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, 0.32, 'lose switch', 'horizontalalignment', 'center', 'fontsize', fz);

text(0, -0.32, 'win switch', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, -0.36, 'lose stay', 'horizontalalignment', 'center', 'fontsize', fz);

text(0.36, .05, 'stay', 'rotation', 270, 'fontsize', fz);
text(-0.36, -.06, 'switch', 'rotation', 90, 'fontsize', fz);

% layout
maxlim = 0.5;
xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
maxlim = 0.4;
set(gca, 'xtick', -maxlim:maxlim:maxlim, 'ytick', -maxlim:maxlim:maxlim);
xlabel('Choice weight'); ylabel('Stimulus weight');
box on; axis square;

end