function [] = fig4g_SjCorrelation(lagGroups, whichmodulator)
% correct - error = stim + resp --stim +- resp = 2stim
% error - correct = -stim + resp -resp - stim = -2stim

% correct + error =

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

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
hold on;
bwMat = cat(3, [groupedDat.response_pupil groupedDat.correct_pupil groupedDat.incorrect_pupil]);
differenceScore = bwMat(:, 2) - bwMat(:, 3); % = 2 * stimweight;

% % write to table
% mat = [bwMat(:, 2); bwMat(:, 3)];
% sjs = dat.response(:, 1) > 0;
% correct = [ones(length(sjs), 1); zeros(length(sjs), 1)];
% tb = array2table([[1:27 1:27]' mat correct [sjs; sjs]] , 'VariableNames', {'subject', 'beta', 'correctness', 'sjGroup'});
% writetable(tb, sprintf('~/Data/pupilUncertainty/CSV/sequentialEffects.csv'));

for sp = 2,
    
    %differenceScore = dat.stimulus_pupil(:, 1); % = 2 * stimweight;
    % differenceScore = bwMat(:, sp); % = 2 * stimweight;
    
    % split subgroups by plain weights
    load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
    
    repSj = find(dat.response(:, 1) > 0);
    switchSj = find(dat.response(:, 1) < 0);    
    % subplot(4,4,sp-1); hold on;
    
    scatter(dat.response(:, 1), differenceScore, 1, 'k', 'filled');
    l = lsline; l.Color = colors(9, :);
    
    scatter(dat.response(repSj, 1), differenceScore(repSj), 10, colors(2, :), 'filled');
    scatter(dat.response(switchSj, 1), differenceScore(switchSj), 10, colors(5, :), 'filled');
    
    [rho, pval] = corr(dat.response(:, 1), differenceScore, 'type', 'spearman');
    text(0.15, -0.15, sprintf('r = %.3f', rho));
    text(0.15, -0.2, sprintf('p < %.3f', pval));
    ylim([-.25 .25]); set(gca, 'ytick', [-.2 0 0.2]);
    xlim([-0.5 0.5]); set(gca, 'xtick', [-.4 0 0.4]);
    xlabel('Response weight');
    
    if sp == 2,
        ylabel({'Pupil modulation'; 'of correct weights'});
    else
        ylabel({'Pupil modulation'; 'of error weights'});
    end
    
    ylabel({'Delta pupil modulation'; 'correct vs. error weights'});
    axis square; box off;
end



%% test a 3-way correlation
cmap = cbrewer('div', 'PuOr', 27);
colormap(cmap);
[val, idx] = sort(dat.response(:, 1));
scatter3(dat.response(:, 1), bwMat(:, 2), bwMat(:, 3), 'filled');
xlabel('Response weight'); ylabel('Correct'); zlabel('Error');