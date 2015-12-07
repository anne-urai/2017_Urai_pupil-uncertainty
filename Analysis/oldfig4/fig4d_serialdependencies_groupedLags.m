% ============================================ %
% correlation between history effect and pupil interaction
% ============================================ %
% group lags 4-7

clear; clc;
whichmodulator = 'pupil';
subjects = 1:27;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

% how do we want to group the lags?
lagGroups = {[1], [1:2], [1:3], [4:7]};

names = fieldnames(dat);
names = names(1:end-2); % remove the two pvals
for i = 1:length(names),
    olddat = dat.(names{i});
    newdat(:, 1:3) = olddat(:, 1:3);
    newdat(:, 4) = mean(olddat(:, 4:end), 2);
    % dont mean 4-7 here!
    datG.(names{i}) = olddat;
end

for l = 1:length(lagGroups)
    subplot(4,4,l); hold on;
    if l ==  1, ylabel('Response weight'); end
    
      % start with a graph for all the trials combined
    bar(1, squeeze(mean(mean(datG.response(:, lagGroups{l}), 2))), ...
        'facecolor', [0.6 0.6 0.6], 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(1, squeeze(mean(mean(datG.response(:, lagGroups{l}), 2))), ...
        [], squeeze(std(mean(datG.response(:, lagGroups{l}), 2))) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    axis tight; axis square;
    
    [pval] = randtest1d(mean(datG.response(:, lagGroups{l}), 2));
    sigstar({[1 1]}, pval);
    title(sprintf('p = %.3f', pval));

    xlabel(sprintf('Lags %d, %d, %d, %d', lagGroups{l}));
    
end

% also for pupil
for l = 1:length(lagGroups)
    subplot(4,4,l+8); hold on;
    if l ==  1, ylabel(sprintf('Response*%s weight', whichmodulator)); end
    
       % start with a graph for all the trials combined
    bar(1, squeeze(mean(mean(datG.response_pupil(:, lagGroups{l}), 2))), ...
        'facecolor', [0.6 0.6 0.6], 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(1, squeeze(mean(mean(datG.response_pupil(:, lagGroups{l}), 2))), ...
        [], squeeze(std(mean(datG.response_pupil(:, lagGroups{l}), 2))) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
   % axis tight; 
    axis square; box off;
    
    [pval] = randtest1d(mean(datG.response_pupil(:, lagGroups{l}), 2));
    sigstar({[1 1]}, pval);

    xlabel(sprintf('Lags %d, %d, %d, %d', lagGroups{l}));
    title(sprintf('p = %.3f', pval));
     
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/history_groupedLags_%s.pdf', whichmodulator));
