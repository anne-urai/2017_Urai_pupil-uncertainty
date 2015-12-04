% ============================================ %
% correlation between history effect and pupil interaction
% ============================================ %
% group lags 4-7

clear; clc;
whichmodulator = 'pupil';
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

% how do we want to group the lags?
lagGroups = {[1], [2:3], [4:7]};

names = fieldnames(dat);
names = names(1:end-2); % remove the two pvals
for i = 1:length(names),
    olddat = dat.(names{i});
    newdat(:, 1:3) = olddat(:, 1:3);
    newdat(:, 4) = mean(olddat(:, 4:end), 2);
    datG.(names{i}) = olddat;
end

for l = 1:length(lagGroups)
    subplot(4,4,l); hold on;
    if l ==  1, ylabel('Stimulus weight'); end
    
    plot(mean(datG.response(:, lagGroups{l}), 2), mean(datG.stimulus(:, lagGroups{l}), 2), ...
        'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
    axis tight; axis square; xlims = get(gca, 'xlim'); ylims = get(gca, 'ylim');
    
    plot([-1 1], [-1 1], 'color', [0.5 0.5 0.5]);
    plot([-1 1], [1 -1], 'color', [0.5 0.5 0.5]);
    
    maxlim = 1.1*max(abs([xlims ylims]));
    maxlim = 0.3;
    xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
    title(sprintf('Lags %d, %d, %d, %d', lagGroups{l}));
    
end
suplabel('Response weight', 'x');

% also for pupil
for l = 1:length(lagGroups)
    subplot(4,4,l+8); hold on;
    if l ==  1, ylabel(sprintf('Stimulus*%s weight', whichmodulator)); end
    
    plot(mean(datG.response_pupil(:, lagGroups{l}), 2), mean(datG.stimulus_pupil(:, lagGroups{l}), 2), ...
        'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
    axis tight; axis square; xlims = get(gca, 'xlim'); ylims = get(gca, 'ylim');
    
    plot([-1 1], [-1 1], 'color', [0.5 0.5 0.5]);
    plot([-1 1], [1 -1], 'color', [0.5 0.5 0.5]);
    
    maxlim = 1.1*max(abs([xlims ylims]));
    maxlim = 0.15;
    xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
    title(sprintf('Lags %d, %d, %d, %d', lagGroups{l}));
    
end
suplabel(sprintf('Response*%s weight', whichmodulator), 'x');
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/history_strategies_%s.pdf', whichmodulator));
