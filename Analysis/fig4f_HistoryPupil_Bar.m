function fig4f_HistoryPupil_Bar(lagGroups)

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

clc; 
whichmodulator = 'pupil';
subjects = 1:27;
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

subplot(4,4,13); hold on;
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
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4e_historyPupil_%s.pdf', whichmodulator));
