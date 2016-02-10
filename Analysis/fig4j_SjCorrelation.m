function [] = fig4j_SjCorrelation(whichmodulator)
% plot correlation between subjects

global mypath;
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));



% or without errorbars
scatter(dat.response(:, 1), dat.(['response_' whichmodulator])(:, 1), 10, mycolmap, 'filled');

[rho, pval] = corr(dat.response(:, 1), dat.(['response_' whichmodulator])(:, 1), 'type', 'pearson');
max(abs(dat.response(:, 1)))
max(abs(dat.(['response_' whichmodulator])(:, 1)))
if pval < 0.05,
   lsline;
end
axis square; box on;

% plot with errorbars
for sj = 1:27, 
    hold on;
    h = ploterr(dat.response(sj, 1), dat.(['response_' whichmodulator])(sj, 1), ...
        {dat.responseCI(sj, 1, 1) dat.responseCI(sj, 1, 2)}, ...
        {dat.(['response_' whichmodulator 'CI'])(sj, 1, 1) ... 
        dat.(['response_' whichmodulator 'CI'])(sj, 1, 2)}, '.', 'abshhxy', 0);
    set(h(1), 'color', mycolmap(sj, :), 'markerfacecolor', mycolmap(sj, :), 'markersize', 1);
    set(h(2), 'color', mycolmap(sj, :), 'linewidth', 0.5);
    set(h(3), 'color', mycolmap(sj, :), 'linewidth', 0.5);
end

text(0.5, -0.15, sprintf('r = %.3f', rho));
text(0.5, -0.2, sprintf('p < %.3f', pval));
ylim([-.25 .25]); set(gca, 'ytick', [-.2 0 0.2]);
xlim([-0.5 0.5]); set(gca, 'xtick', [-.4 0 0.4]);

bf10 = corrbf(rho,27);
text(0.5, -0.25, sprintf('bf10 = %.3f', bf10));

switch whichmodulator
    case 'pupil'
        ylabel('Pupil');
    case 'rt'
        ylabel('RT');
        xlabel('Choice weight');
        
        % between subject stuff
        x = dat.response(:, 1);
        y1 = dat.response_pupil(:, 1);
        y2 = dat.response_rt(:, 1);
        [pval, deltaR] = permtest_correlation(x, y1, y2, 0, 1000);
end
set(gca, 'xcolor', 'k', 'ycolor', 'k');

end
