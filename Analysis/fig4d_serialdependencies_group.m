function [pvalS, pvalR] = fig4d_serialdependencies_group

clear all; close all; clc;

whichmodulator = 'rt';
load(sprintf('~/Data/UvA_pupil/GrandAverage/historyweights_%s.mat', whichmodulator));

% ============================================ %
% grand average history kernels
% ============================================ %

clf;
lags = 1:7; subjects = 27;
colors = linspecer(8);
nlags = 7;
ylim1st = [-0.4 0.4];
ylim2nd = [-0.15 0.15];

clf;
% stim resp
subplot(3,3,1);
bh = boundedline(lags, nanmean(dat.stimulus), nanstd(dat.stimulus) ./ sqrt(length(subjects)), ...
    lags, nanmean(dat.response), nanstd(dat.response) ./ sqrt(length(subjects)), ...
    'cmap', colors([2 4], :), 'alpha');
doSigstar(dat.stimulus, 2); doSigstar(dat.response, 4);
legend(bh, 'stimulus', 'response'); legend boxoff;
xlim([0.5 nlags+0.5]); ylim(ylim1st); set(gca, 'xtick', lags, 'ytick', -1:0.2:1);
ylabel('History kernels');

% err correct
subplot(3,3,2);
bh = boundedline(lags, nanmean(dat.correct), nanstd(dat.correct) ./ sqrt(length(subjects)), ...
    lags, nanmean(dat.incorrect), nanstd(dat.incorrect) ./ sqrt(length(subjects)), ...
    'cmap', colors([3 1], :), 'alpha');
doSigstar(dat.correct, 3); doSigstar(dat.incorrect, 1);
legend(bh, 'correct', 'error'); legend boxoff;
axis tight; ylims = get(gca, 'ylim'); set(gca, 'ylim', [-max(abs(ylims)*1.1) max(abs(ylims)*1.1)]);
xlim([0.5 nlags+0.5]); ylim(ylim1st);set(gca, 'xtick', lags, 'ytick', -1:0.2:1);

% now the pupil
subplot(3,3,4);
bh = boundedline(lags, nanmean(dat.stimulus_pupil), nanstd(dat.stimulus_pupil) ./ sqrt(length(subjects)), ...
    lags, nanmean(dat.response_pupil), nanstd(dat.response_pupil) ./ sqrt(length(subjects)), ...
    'cmap', colors([2 4], :), 'alpha');
pvalS = doSigstar(dat.stimulus_pupil, 2); pvalR = doSigstar(dat.response_pupil, 4);
legend(bh, 'stimulus', 'response'); legend boxoff;
axis tight; ylims = get(gca, 'ylim'); set(gca, 'ylim', [-max(abs(ylims)*1.1) max(abs(ylims)*1.1)]);
xlim([0.5 nlags+0.5]); ylim(ylim2nd);
ylabel([whichmodulator ' interaction']);set(gca, 'xtick', lags, 'ytick', -1:0.1:1);

subplot(3,3,5);
bh = boundedline(lags, nanmean(dat.correct_pupil), nanstd(dat.correct_pupil) ./ sqrt(length(subjects)), ...
    lags, nanmean(dat.incorrect_pupil), nanstd(dat.incorrect_pupil) ./ sqrt(length(subjects)), ...
    'cmap', colors([3 1], :), 'alpha');
doSigstar(dat.correct_pupil, 3); doSigstar(dat.incorrect_pupil, 1);
legend(bh, 'correct', 'error'); legend boxoff;
axis tight; ylims = get(gca, 'ylim'); set(gca, 'ylim', [-max(abs(ylims)*1.1) max(abs(ylims)*1.1)]);
xlim([0.5 nlags+0.5]); ylim(ylim2nd);set(gca, 'xtick', lags, 'ytick', -1:0.1:1);

subplot(3,3,6);
bh = boundedline(lags, nanmean(dat.pupil), nanstd(dat.pupil) ./ sqrt(length(subjects)), ...
    'cmap', colors([5], :), 'alpha');
doSigstar(dat.pupil, 5); 
legend(bh, whichmodulator); legend boxoff;
axis tight; ylims = get(gca, 'ylim'); set(gca, 'ylim', [-max(abs(ylims)*1.1) max(abs(ylims)*1.1)]);
xlim([0.5 nlags+0.5]); ylim(ylim2nd);set(gca, 'xtick', lags, 'ytick', -1:0.1:1);
xlabel('Lags');

suplabel('Grand Average', 't');
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/serial/historykernels_%s_GA.pdf', whichmodulator));

end

function pval = doSigstar(mat, colIdx)

nanidx = find(isnan(nanmean(mat, 2)));
mat(nanidx, :) = [];

colors = linspecer(8);
% for each lag, do permutation test
for m = 1:7,
   [~, pval(m)] = permtest(mat(:, m));
end

signific = find(pval < 0.05);
meanMat = nanmean(mat);
plot(signific, meanMat(signific), ...
    'o', 'markersize', 4, 'markeredgecolor', 'w', 'markerfacecolor', colors(colIdx, :));

end