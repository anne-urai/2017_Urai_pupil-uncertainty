% plot the history kernels, only mean +- sem
global mypath
close all

nbins = 3;
grandavg = postPupilBehaviour('baseline_pupil', nbins, []);
subplot(4,4,1);
% line to indicate 50% repetition
y = grandavg.repetition;
plot([1 nbins], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;
thiscolor = [0 0 0]; thismarker = '.'; thismarkersize = 14;
% errorbar
h = ploterr(1:nbins, nanmean(y), [], nanstd(y) ./sqrt(27), 'k-',  'abshhxy', 0);
set(h(1), 'color', thiscolor, 'markersize', thismarkersize, 'marker',thismarker);
set(h(2), 'color', thiscolor); % line color

xticklabs       = repmat({' '}, 1, nbins);
xticklabs{1}    = 'low';
xticklabs{end}  = 'high';
if nbins == 3, xticklabs{2} = 'med'; end

set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
    'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
% axis square;
xlim([0.5 nbins+0.5]);

% determine y label and limits
set(gca, 'ylim', [0.46 0.56], 'ytick', 0.48:0.02:0.56);
ylabel('P(repeat)');

% do Bayesian ANOVA to get Bayes Factors
statdat             = table;
statdat.DV          = y(:);
sj = repmat(1:27, nbins, 1)';
statdat.subjnr      = sj(:);
ft      = repmat(transpose(1:nbins), 1, 27)';
statdat.prevPupilBins = ft(:);
writetable(statdat, sprintf('%s/Data/CSV/ANOVAdat.csv', mypath));
system('/Library/Frameworks/R.framework/Resources/bin/R < BayesFactorANOVA.R --no-save');
statres = readtable(sprintf('%s/Data/CSV/ANOVAresults.csv', mypath)); % fetch results

yval = max(get(gca, 'ylim'));
if statres.pvalue < 0.05, % only show if significant
    mysigstar(gca, [1 nbins], [yval yval], statres{cnt}.pvalue, 'k', 'down');
end
xlabel('Baseline pupil'); axis square;

subplot(442);
plot([1 7], [0 0], ':'); hold on;
set(gca, 'ycolor', 'k', 'xtick', 1:7, 'ytick', [-0.1:0.05:0.1], ...
    'ylim', [-.07 .01], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
ylabel('Pupil x choice weight');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
xlabel('Lags'); axis square;
plot(1, -0.06, '*k', 'markersize', 2);
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
boundedline(1:7, mean(dat.response_pupil), std(dat.response_pupil) ./ sqrt(27), 'k');

print(gcf, '-dpdf', sprintf('%s/Figures/S7_pupilLags.pdf', mypath));
