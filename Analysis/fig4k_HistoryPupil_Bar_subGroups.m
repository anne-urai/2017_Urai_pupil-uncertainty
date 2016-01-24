function fig4k_HistoryPupil_Bar_subGroups(whichmodulator, grouping)

if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('grouping', 'var'); grouping = 'repeat'; end

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'pupil-rt'));
pupMat = dat.response_pupil(:, 1);

load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'rt-pupil'));
rtMat = dat.response_pupil(:, 1);

% ============================================ %
% barweb matrix
% ============================================ %

colors = cbrewer('qual', 'Set1', 9);
bwMat = cat(3, pupMat, rtMat);

% split subgroups by plain weights
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));

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
[~, pvals(3)] = ttest(bwMat(theseSj, 1), bwMat(theseSj, 2));

hold on;
barcolors = colors([9 9], :);

for i = 1:2,
    bar(i, mean(bwMat(theseSj, i)), 'barwidth', 0.5', 'facecolor', barcolors(i, :), 'edgecolor', 'none');
    errorbar(i, mean(bwMat(theseSj, i)), std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)), 'k');
    % add significance star
    
    if mean(bwMat(theseSj, i)) < 0,
        yval = min(bwMat(theseSj, i)) - 1*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)));
    else
        yval = mean(bwMat(theseSj, i)) + 2*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)));
    end
    yval = -0.12;
    mysigstar(i, yval, pvals(i), barcolors(i, :));
end

mysigstar([1 2], 0.01, pvals(3), 'k');

ylims = get(gca, 'ylim');
set(gca, 'ylim', [ylims(1) - 0.2*range(ylims) ylims(2)+0.2*range(ylims)]);
ylim([-0.12 0.02]); set(gca, 'ytick', [-0.1:0.1:0]);
set(gca, 'xtick', 1:2, ...
    'xticklabel', {'pupil', 'rt'}, 'xticklabelrotation', 0, ...
    'xaxislocation', 'top');
ylabel({'Modulation weights'});
axis square;
titlecols = cbrewer('div', 'PuOr', 10);

switch grouping
    case 'switch'
        title('Alternators', 'color', titlecols(1, :));
    case 'repeat'
        title('Repeaters', 'color', titlecols(end, :));
end

% write to table for JASP
mat = [pupMat rtMat];
sjs = dat.response(:, 1) > 0;
tb = array2table([transpose(1:27) mat sjs] , 'VariableNames', {'sjnr', 'pupil', 'rt', 'sjGroup'});
writetable(tb, sprintf('~/Data/pupilUncertainty/CSV/sequentialEffects.csv'));

end
