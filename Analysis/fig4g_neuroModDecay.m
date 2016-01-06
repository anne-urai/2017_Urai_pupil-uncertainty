function [] = fig4g_neuroModDecay(lagGroups, whichmodulator, grouping, correctness)

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'baseline'; end
if ~exist('grouping', 'var'); grouping = 'all'; end
if ~exist('correctness', 'var'); correctness = 1; end

% ========================================================= %
% panel E: neuromodulatory decay
% ========================================================= %

nlags = 30;
hold on;

subjects = 1:27;
for sj = unique(subjects),
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    %  trls = 1:size(data, 1);
    trls = find(data.correct == correctness);
    
    
    % first, get the correlation at lag 0 between decision pupil and the
    % end of the trial
    tmp = corrcoef(data.decision_pupil(trls)', data.trialend_pupil(trls)', 'rows', 'complete');
    neuromodDecay.sametrial(sj) = tmp(2);
    
    for lag = 1:nlags,
        
        % with this selection, take the trials after that
        laggedtrls = trls+lag;
        
        % exclude trials at the end of the block
        if any(laggedtrls > size(data, 1)),
            trls(laggedtrls > size(data, 1)) = [];
            laggedtrls(laggedtrls > size(data, 1)) = [];
        end
        
        % remove trials that dont match in block nr
        noblockmatch = data.blocknr(laggedtrls) ~= data.blocknr(trls);
        laggedtrls(noblockmatch) = [];
        trls(noblockmatch) = [];

        % these look more or less the same
        neuromodDecay.pear(sj, lag) = corr(data.decision_pupil(trls), data.baseline_pupil(laggedtrls));
    end
end


% split subjects based on their plain history weights
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
colors = cbrewer('qual', 'Set1', 9);
switch grouping
    case 'all'
        theseSj = 1:27;
        tit = 'All subjects';
        titcolor = 'k';
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
        titcolor = colors(2,:);
        tit = 'Repeaters';
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
        titcolor = colors(5,:);
        tit = 'Switchers';
end

colors = linspecer(3);
switch correctness
    case 1
        colors = colors(3, :);
    case 0
        colors = colors(2, :);
end

hold on;
plot([1 lag], [0 0], 'k', 'linewidth', 0.5);

errorbar(0, nanmean(neuromodDecay.sametrial(theseSj)), (nanstd(neuromodDecay.sametrial(theseSj)) ./ sqrt(length(subjects))), ...
    '.', 'Color', colors, 'MarkerFaceColor', colors, 'markersize', 8);

boundedline(1:lag, mean(neuromodDecay.pear(theseSj, :, 1)), std(neuromodDecay.pear(theseSj, :, 1)) ./ sqrt(length(subjects)), 'cmap', colors);
axis tight;  xlim([0 lag+0.5]);
xlabel('Lags'); ylabel({'Pearson''s rho'});
set(gca, 'xtick', [1 nlags/2 nlags]);
ylim([-0.08 0.12]);
set(gca, 'xcolor', 'k', 'ycolor', 'k');
axis square;
title(tit, 'color', titcolor);
