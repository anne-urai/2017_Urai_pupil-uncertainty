function fig3hi_HistoryPupil_Bar(whichmodulator)
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
idx = 1:2;

% p values
[~, pvals(1)] = ttest(bwMat(theseSj, 1));
[~, pvals(2)] = ttest(bwMat(theseSj, 2));
[~, pvals(3)] = ttest(bwMat(theseSj, 1), bwMat(theseSj, 2));
disp(pvals);

hold on;
barcolors = colors([9 9 3 1], :);

for i = idx,
    bar(find(idx==i), mean(bwMat(theseSj, i)), 'barwidth', 0.5', 'facecolor', barcolors(i, :), 'edgecolor', 'none');
    h = ploterr(find(idx==i), mean(bwMat(theseSj, i)), [], ...
        std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)), 'k', 'abshhxy', 0);
    set(h(1), 'marker', 'none');
    % add significance star
    if mean(bwMat(theseSj, i)) < 0,
        yval = min(bwMat(theseSj, i)) - 1*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)));
    else
        yval = mean(bwMat(theseSj, i)) + 2*(std(bwMat(theseSj, i)) ./ sqrt(length(theseSj)));
    end
    yval = -0.1;
    mysigstar(find(idx==i), yval, pvals(i));
end

yval = -0.12;
% sigstar between the two
mysigstar([1 2], yval, pvals(3));

ylims = get(gca, 'ylim');
set(gca, 'ylim', [ylims(1) - 0.2*range(ylims) ylims(2)+0.2*range(ylims)]);
ylim([-0.12 0.02]); set(gca, 'ytick', [-0.1:0.1:0]);
set(gca, 'xtick', 1:4, 'xticklabel', [], 'xaxislocation', 'top');
ylabel({whichmodulator; 'modulation weights'});
% axis square;

end
