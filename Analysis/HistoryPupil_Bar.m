function HistoryPupil_Bar(whichmodulator)
global mypath;

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

switch whichmodulator
    case 'pupil'
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
    case 'rt'
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
    otherwise
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator));
        whichmodulator = 'pupil'; % use the first modulatory term
end

% ============================================ %
% barweb matrix
% ============================================ %

colors = cbrewer('qual', 'Set1', 9);
bwMat = cat(3, [dat.(['response_' whichmodulator])(:, 1) dat.(['stimulus_' whichmodulator])(:, 1)]);

% split subgroups by plain weights
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));

theseSj = 1:27;
% p values
[~, pvals(1)] = ttest(bwMat(theseSj, 1));
[~, pvals(2)] = ttest(bwMat(theseSj, 2));
[~, pvals(3)] = ttest(bwMat(theseSj, 1), bwMat(theseSj, 2));

hold on;
barcolors = colors([9], :);
i = [1 2];
bar(i, mean(bwMat(theseSj, i)), 'barwidth', 0.5', 'facecolor', barcolors, 'edgecolor', 'none');
h = ploterr(i, mean(bwMat(theseSj, i)), [], ...
    std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)), 'k', 'abshhxy', 0);
set(h(1), 'marker', 'none');
plot(bwMat(theseSj, i)', '.k-', 'linewidth', 0.2);

% sigstar between the two
mysigstar(gca, 1, 0.2, pvals(1));
mysigstar(gca, 2, 0.2, pvals(2));
mysigstar(gca, [1 2], -0.4, pvals(3));

% ylims = get(gca, 'ylim');
% set(gca, 'ylim', [ylims(1) - 0.2*range(ylims) ylims(2)+0.2*range(ylims)]);
ylim([-0.4 0.2]); set(gca, 'ytick', [-0.4:0.2:0.2]);
set(gca, 'xtick', i, 'xticklabel', [], 'xaxislocation', 'top');
xlim([0.5 2.5]);
set(gca, 'xcolor', 'w', 'ycolor', 'k');

% already make the titles right
titcols = cbrewer('div', 'RdYlBu', 12);
switch whichmodulator
    case 'pupil'
        ylabel({'Pupil response'; 'modulation weights'}, 'color', titcols(1,:));
    case 'rt'
        ylabel({'Reaction time'; 'modulation weights'}, 'color', titcols(end, :));
        set(gca, 'yaxislocation', 'right');
        set(get(gca, 'ylabel'), 'rotation', 270);
end

end
