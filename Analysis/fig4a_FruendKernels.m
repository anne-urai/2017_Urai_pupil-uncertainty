function [] = fig4a_FruendKernels(whichmodulator)

% make all the panels for figure 4a
%     - panel A: frund kernels
%     - panel B: decision strategy for lags 1-3
%     - panel C: valueShift on lag 1-3 with complete stats anova
%     - panel D: frund interaction lags 1-3, explain why we use the multiplicative term
%     - panel D: lag 1 valueShift + Frund interaction bargraph
%     - panel E: decay of neuromodulation over trials (baseline correct the next 7 baselines)
%     - panel F: correct vs error, model comparison

addpath('~/Documents/fieldtrip');
ft_defaults;

% determine the subjects based on their plain weights
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

% save colors
colors = cbrewer('div', 'PuOr', 256);
%colors(1:128, :) = flipud(colors(129:end, :));
% flip the halves
colors(128-40:128+40, :) = []; % remove white in the middle
colors = colors + 0.1; % make a bit paler
colors(colors > 1) = 1;

for sj = 1:27,
    % find which color this is
    colspace =  linspace(-0.5, 0.5, size(colors, 1));
    colidx = dsearchn(colspace', dat.response(sj, 1));
    mycolmap(sj, :) = colors(colidx, :);
end
save('~/Data/pupilUncertainty/GrandAverage/sjcolormap.mat', 'mycolmap');

hold on;
for sj = 1:27,
    plot(dat.response(sj, :)', 'color', mycolmap(sj, :), 'linewidth', 0.5);
end
scatter(ones(1, 27), dat.response(:, 1), 10, mycolmap, 'filled');

% add the group
[ax, p1, p2] = plotyy(1:7, nanmean(dat.response),  ...
    1:7, mean(abs(dat.response)));

set(ax(1), 'ycolor', 'k', 'xtick', 1:7, 'ytick', [-0.4 0 0.4], 'ylim', [-.45 .4], 'xlim', [0.5 7.5], 'box', 'off');
set(ax(2), 'ycolor', [0.4 0.4 0.4], 'xtick', 1:7, 'ytick', [0 0.2], 'ylim', [-.2 .2], 'xlim', [0.5 7.5], 'box', 'off');

set(p1, 'color', 'k', 'linewidth', 1);
set(p2, 'color', [0.4 0.4 0.4], 'linewidth', 1);
axis(ax(1), 'square'); axis(ax(2), 'square');

ylabel(ax(1),'Choice weight') % label left y-axis
ylabel(ax(2),'|Choice weight|') % label right y-axis
xlabel('Lags');
end
