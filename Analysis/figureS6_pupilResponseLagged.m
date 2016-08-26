% plot the history kernels, only mean +- sem
global mypath
close all
subplot(441);
plot([1 7], [0 0], ':'); hold on;
set(gca, 'ycolor', 'k', 'xtick', 1:7, 'ytick', [-0.1:0.05:0.1], ...
    'ylim', [-.07 .01], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
ylabel('Pupil x choice weight');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
xlabel('Lags'); axis square;
plot(1, -0.06, '*k', 'markersize', 2);
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
boundedline(1:7, mean(dat.response_pupil), std(dat.response_pupil) ./ sqrt(27), 'k');

print(gcf, '-dpdf', sprintf('%s/Figures/pupilLags.pdf', mypath));
