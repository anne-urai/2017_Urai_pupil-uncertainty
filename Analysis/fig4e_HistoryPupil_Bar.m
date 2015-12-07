% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

clear; clc; clf;
whichmodulator = 'pupil';
subjects = 1:27;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

subplot(4,4,1); hold on;
% start with a graph for all the trials combined
bar(1, squeeze(mean(dat.response(:, 1))), ...
    'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none', 'barwidth', 0.5);
h = ploterr(1, squeeze(mean(dat.response(:, 1))), ...
    [], squeeze(std(dat.response(:, 1))) / sqrt(length(subjects)), ...
    'k.', 'hhxy', 0.0000000000000001);
set(h(1), 'marker', 'none');

% also pupil
bar(2, squeeze(mean(dat.response_pupil(:, 1))), ...
    'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none', 'barwidth', 0.5);
h = ploterr(2, squeeze(mean(dat.response_pupil(:, 1))), ...
    [], squeeze(std(dat.response_pupil(:, 1))) / sqrt(length(subjects)), ...
    'k.', 'hhxy', 0.0000000000000001);
set(h(1), 'marker', 'none');

% add stats
[pval(1)] = signrank(dat.response(:, 1));
[pval(2)] = signrank(dat.response_pupil(:, 1));

sigstar({[1 1], [2 2]}, pval);
set(gca, 'xtick', 1:2, 'xticklabel', {'resp_{-1}', sprintf('resp_{-1}*%s_{-1}', whichmodulator)}, 'xticklabelrotation', -20);
ylabel('Beta weights');


%% now correct trials
colors = linspecer(3); colors = colors(3, :);
subplot(4,4,2); hold on;
% start with a graph for all the trials combined
bar(1, squeeze(mean(dat.correct(:, 1))), ...
    'facecolor', colors, 'edgecolor', 'none', 'barwidth', 0.5);
h = ploterr(1, squeeze(mean(dat.correct(:, 1))), ...
    [], squeeze(std(dat.correct(:, 1))) / sqrt(length(subjects)), ...
    'k.', 'hhxy', 0.0000000000000001);
set(h(1), 'marker', 'none');

% also pupil
bar(2, squeeze(mean(dat.correct_pupil(:, 1))), ...
    'facecolor', colors, 'edgecolor', 'none', 'barwidth', 0.5);
h = ploterr(2, squeeze(mean(dat.correct_pupil(:, 1))), ...
    [], squeeze(std(dat.correct_pupil(:, 1))) / sqrt(length(subjects)), ...
    'k.', 'hhxy', 0.0000000000000001);
set(h(1), 'marker', 'none');

% add stats
[pval(1)] = signrank(dat.correct(:, 1));
[pval(2)] = signrank(dat.correct_pupil(:, 1));

sigstar({[1 1], [2 2]}, pval);
set(gca, 'xtick', 1:2, 'xticklabel', {'resp_{-1}', sprintf('resp_{-1}*%s_{-1}', whichmodulator)}, 'xticklabelrotation', -20);
title('Correct trials');

ylabel('Beta weights');
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4e_historyPupil_%s.pdf', whichmodulator));
