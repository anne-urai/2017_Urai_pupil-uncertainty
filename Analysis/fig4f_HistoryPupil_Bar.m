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

colors = linspecer(9);
bwMat = cat(3, [groupedDat.response_pupil groupedDat.correct_pupil groupedDat.incorrect_pupil]);

switch grouping
    
    case 'all'
        theseSj = 1:27;
        tit = 'All subjects';
        thiscolor = colors(2, :);
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
        thiscolor = colors(8, :);
        tit = 'Repeaters';
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
        thiscolor = colors(9, :);
        tit = 'Switchers';
end

% sigstars
pvals = squeeze(cat(3, ...
    [randtest1d(groupedDat.response_pupil(theseSj)) randtest1d(groupedDat.correct_pupil(theseSj)) ...
    randtest1d(groupedDat.incorrect_pupil(theseSj))]));
hold on;

barcolors = colors([2 3 1], :);

for i = 1:3,
    bar(i, mean(bwMat(theseSj, i)), 'barwidth', 0.5', 'facecolor', barcolors(i, :), 'edgecolor', 'none');
    errorbar(i, mean(bwMat(theseSj, i)), std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)), 'k');
    
    % add significance star
    if mean(bwMat(theseSj, i)) < 0,
        mysigstar(i, mean(bwMat(theseSj, i)) - 2*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj))), pvals(i));
    else
        mysigstar(i, mean(bwMat(theseSj, i)) + 2*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj))), pvals(i));
    end
    
end

switch grouping
    case 'all'
        ylabel('Resp-1*Pupil-1');
end
set(gca, 'xtick', []);
ylims = get(gca, 'ylim');
set(gca, 'ylim', [ylims(1) - 0.2*range(ylims) ylims(2)+0.2*range(ylims)]);
%axis tight; axis square; xlim([0.5 i+0.5]);
title(tit, 'color', thiscolor);

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4e_historyPupil_%s.pdf', whichmodulator));
