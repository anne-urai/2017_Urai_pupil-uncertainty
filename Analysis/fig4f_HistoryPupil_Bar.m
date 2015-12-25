function fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, grouping)

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('grouping', 'var'); grouping = 'all'; end

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

clc;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

% ============================================ %
% combine lags
% ============================================ %
flds = fieldnames(dat);

for f = 1:length(flds),
    try
        groupedDat.(flds{f}) = mean(dat.(flds{f})(:, lagGroups), 2);
    end
end

% ============================================ %
% barweb matrix
% ============================================ %

colors = cbrewer('qual', 'Set1', 9);
bwMat = cat(3, [groupedDat.response_pupil groupedDat.correct_pupil groupedDat.incorrect_pupil]);

% split subgroups by plain weights
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));

switch grouping
    
    case 'all'
        theseSj = 1:27;
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
end

% sigstars
pvals = squeeze(cat(3, ...
    [permtest(groupedDat.response_pupil(theseSj)) permtest(groupedDat.correct_pupil(theseSj)) ...
    permtest(groupedDat.incorrect_pupil(theseSj))]));
hold on;

barcolors = colors([9 3 1], :);

for i = 1:3,
    bar(i, mean(bwMat(theseSj, i)), 'barwidth', 0.5', 'facecolor', barcolors(i, :), 'edgecolor', 'none');
    errorbar(i, mean(bwMat(theseSj, i)), std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)), 'k');
    
    % add significance star
    if mean(bwMat(theseSj, i)) < 0,
        mysigstar(i, mean(bwMat(theseSj, i)) - 2*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj))), pvals(i), barcolors(i, :));
    else
        mysigstar(i, mean(bwMat(theseSj, i)) + 2*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj))), pvals(i));
    end
    
end

% difference error and correct
[pval] = permtest(bwMat(theseSj, 2), bwMat(theseSj, 3));
ymax = min( nanmean(bwMat(theseSj, 2:3)) - ...
    3* nanstd(bwMat(theseSj, 2:3)) ./ sqrt(length(theseSj)));
mysigstar([2 3], [ymax ymax], pval, 'k', 'up');

switch grouping
    case 'all'
        ylabel({'Pupil modulation'; 'of response weights'});
end

set(gca, 'xtick', []);
ylims = get(gca, 'ylim');
set(gca, 'ylim', [ylims(1) - 0.2*range(ylims) ylims(2)+0.2*range(ylims)]);
%axis tight; axis square; xlim([0.5 i+0.5]);
%title(tit, 'color', thiscolor);
ylim([-0.1 0.02]); set(gca, 'ytick', [-0.1:0.1:0]);
set(gca, 'xcolor', 'w');
axis square;

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4e_historyPupil_%s.pdf', whichmodulator));
