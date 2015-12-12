function fig4f_HistoryPupil_Bar(lagGroups, whichmodulator)

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

clc;
subjects = 1:27;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));
if ~exist('lagGroups', 'var'), lagGroups = 1; end

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

posRespSj = find(mean(dat.response(:, lagGroups), 2) > 0);
negRespSj = find(mean(dat.response(:, lagGroups), 2) < 0);
colors = linspecer(9);

bwMat = cat(3, [groupedDat.response groupedDat.correct groupedDat.incorrect], ...
    [groupedDat.response_pupil groupedDat.correct_pupil groupedDat.incorrect_pupil]);

% sigstars
pvals.posRespSj = squeeze(cat(3, [randtest1d(groupedDat.response(posRespSj)) randtest1d(groupedDat.correct(posRespSj)) randtest1d(groupedDat.incorrect(posRespSj))], ...
    [randtest1d(groupedDat.response_pupil(posRespSj)) randtest1d(groupedDat.correct_pupil(posRespSj)) randtest1d(groupedDat.incorrect_pupil(posRespSj))]))';
pvals.negRespSj = squeeze(cat(3, [randtest1d(groupedDat.response(negRespSj)) randtest1d(groupedDat.correct(negRespSj)) randtest1d(groupedDat.incorrect(negRespSj))], ...
    [randtest1d(groupedDat.response_pupil(negRespSj)) randtest1d(groupedDat.correct_pupil(negRespSj)) randtest1d(groupedDat.incorrect_pupil(negRespSj))]))';

subplot(4,4,10); hold on;
barwebQH(squeeze(nanmean(bwMat(posRespSj, :, :)))', squeeze(nanstd(bwMat(posRespSj, :, :)))' ./sqrt(length(posRespSj)), ...
    [], 0.5, colors([2 3 1], :), pvals.posRespSj);
set(gca, 'xtick', 1:2, 'xticklabel', {'Resp-1', 'Resp-1*Pup-1'});
%axis tight; 
ylim([-.15 0.3]);
axis square; xlim([0.5 2.5]);
title('Repeaters', 'color', colors(8,:));

subplot(4,4,11); hold on;
barwebQH(squeeze(nanmean(bwMat(negRespSj, :, :)))', squeeze(nanstd(bwMat(negRespSj, :, :)))' ./sqrt(length(negRespSj)), ...
    [], 0.5, colors([2 3 1], :), pvals.negRespSj);
title('Switchers', 'color', colors(9,:));
axis square; xlim([0.5 2.5]);
ylim([-.35 0.1]);
set(gca, 'xtick', 1:2, 'xticklabel',{'Resp-1', 'Resp-1*Pup-1'});


if 0,
    % start with a graph for all the trials combined
    bar(1, squeeze(mean(groupedDat.response)), ...
        'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(1, squeeze(mean(groupedDat.response)), ...
        [], squeeze(std(groupedDat.response)) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    
    % also pupil
    bar(2, squeeze(mean(groupedDat.response_pupil)), ...
        'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(2, squeeze(mean(groupedDat.response_pupil)), ...
        [], squeeze(std(groupedDat.response_pupil)) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    
    % add stats
    [pval(1)] = randtest1d(groupedDat.response);
    [pval(2)] = randtest1d(groupedDat.response_pupil);
    
    sigstar({[1 1], [2 2]}, pval);
    set(gca, 'xtick', 1:2, 'xticklabel', {'resp_{-1}', sprintf('resp_{-1}*%s_{-1}', whichmodulator)}, 'xticklabelrotation', -20);
    ylabel('Beta weights');
    xlabel('All trials');
    
    
    %% now correct trials
    colors = linspecer(3); colors = colors(3, :);
    subplot(4,4,14); hold on;
    % start with a graph for all the trials combined
    bar(1, squeeze(mean(groupedDat.correct)), ...
        'facecolor', colors, 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(1, squeeze(mean(groupedDat.correct)), ...
        [], squeeze(std(groupedDat.correct)) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    
    % also pupil
    bar(2, squeeze(mean(groupedDat.correct_pupil)), ...
        'facecolor', colors, 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(2, squeeze(mean(groupedDat.correct_pupil)), ...
        [], squeeze(std(groupedDat.correct_pupil)) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    
    % add stats
    [pval(1)] = randtest1d(groupedDat.correct);
    [pval(2)] = randtest1d(groupedDat.correct_pupil);
    
    sigstar({[1 1], [2 2]}, pval);
    set(gca, 'xtick', 1:2, 'xticklabel', {'resp_{-1}', sprintf('resp_{-1}*%s_{-1}', whichmodulator)}, 'xticklabelrotation', -20);
    xlabel('Correct trials');
    
    ylabel('Beta weights');
    
    %%
    colors = linspecer(9); colors = colors(1, :);
    subplot(4,4,15); hold on;
    % start with a graph for all the trials combined
    bar(1, squeeze(mean(groupedDat.incorrect)), ...
        'facecolor', colors, 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(1, squeeze(mean(groupedDat.incorrect)), ...
        [], squeeze(std(groupedDat.incorrect)) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    
    % also pupil
    bar(2, squeeze(mean(groupedDat.incorrect_pupil)), ...
        'facecolor', colors, 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(2, squeeze(mean(groupedDat.incorrect_pupil)), ...
        [], squeeze(std(groupedDat.incorrect_pupil)) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    
    % add stats
    [pval(1)] = randtest1d(groupedDat.incorrect);
    [pval(2)] = randtest1d(groupedDat.incorrect_pupil);
    
    sigstar({[1 1], [2 2]}, pval);
    set(gca, 'xtick', 1:2, 'xticklabel', {'resp_{-1}', sprintf('resp_{-1}*%s_{-1}', whichmodulator)}, 'xticklabelrotation', -20);
    xlabel('Error trials');
    
    ylabel('Beta weights');
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4e_historyPupil_%s.pdf', whichmodulator));
