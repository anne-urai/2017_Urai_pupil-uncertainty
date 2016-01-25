function fig4hi_HistoryPupil_Bar(whichmodulator, grouping)
if ~exist('grouping', 'var'); grouping = 'all'; end

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator));

% ============================================ %
% barweb matrix
% ============================================ %

colors = cbrewer('qual', 'Set1', 9);
bwMat = cat(3, [dat.response_pupil(:, 1) dat.stimulus_pupil(:, 1) ...
    dat.correct_pupil(:, 1) dat.incorrect_pupil(:, 1)]);

% split subgroups by plain weights
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));

switch grouping
    case 'all'
        theseSj = 1:27;
        repeaters = find(dat.response(:, 1) > 0);
        alternators = find(dat.response(:, 1) < 0);
        idx = 1:4;
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
        idx = [1 3 4];
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
        idx = [1 3 4];
end

% p values
[~, pvals(1)] = ttest(bwMat(theseSj, 1));
[~, pvals(2)] = ttest(bwMat(theseSj, 2));
[~, pvals(3)] = ttest(bwMat(theseSj, 3));
[~, pvals(4)] = ttest(bwMat(theseSj, 4));

hold on;
barcolors = colors([9 9 3 1], :);

for i = idx,
    bar(find(idx==i), mean(bwMat(theseSj, i)), 'barwidth', 0.5', 'facecolor', barcolors(i, :), 'edgecolor', 'none');
    errorbar(find(idx==i), mean(bwMat(theseSj, i)), std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)), 'k');
    % add significance star
    
    if mean(bwMat(theseSj, i)) < 0,
        yval = min(bwMat(theseSj, i)) - 1*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)));
    else
        yval = mean(bwMat(theseSj, i)) + 2*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)));
    end
    yval = -0.12;
    mysigstar(find(idx==i), yval, pvals(i), barcolors(i, :));
end

ylims = get(gca, 'ylim');
set(gca, 'ylim', [ylims(1) - 0.2*range(ylims) ylims(2)+0.2*range(ylims)]);
ylim([-0.12 0.02]); set(gca, 'ytick', [-0.1:0.1:0]);
set(gca, 'xtick', 1:4, ...
    'xticklabel', {'resp', 'stim', 'correct', 'incorrect'}, 'xticklabelrotation', 0, ...
    'xaxislocation', 'top');
ylabel({'Modulation weights'});
axis square;
end
