function [] = SjCorrelation(whichmodulator, whichweight)
% plot correlation between subjects

if ~exist('whichweight', 'var'); whichweight = 'response'; end

global mypath;
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));

% or without errorbars
scatter(dat.response(:, 1), dat.([whichweight '_' whichmodulator])(:, 1), 10, mycolmap, 'filled');

[rho, pval] = corr(dat.response(:, 1), dat.([whichweight '_' whichmodulator])(:, 1), 'type', 'pearson');
if pval < 0.05,
   lsline;
end
axis square; box on;

% plot with errorbars
for sj = 1:27, 
    hold on;
    % the CI fields have absolute bounds, lower and upper
    h = ploterr(dat.response(sj, 1), dat.([whichweight '_' whichmodulator])(sj, 1), ...
        {dat.responseCI(sj, 1, 1) dat.responseCI(sj, 1, 2)}, ...
        {dat.([whichweight '_' whichmodulator 'CI'])(sj, 1, 1) ... 
        dat.([whichweight '_' whichmodulator 'CI'])(sj, 1, 2)}, '.', 'abshhxy', 0);
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
        y1 = dat.([whichweight '_pupil'])(:, 1);
        y2 = dat.([whichweight '_rt'])(:, 1);
        x = dat.(whichweight)(:, 1);
        
        % one way, parametric Steiger's test
        [rddiff,cilohi,p] = rddiffci(corr(x, y1), corr(x, y2), corr(y1, y2), 27, 0.05);
        
        % second way, nonparametric permutation
        [pval, deltaR] = permtest_correlation(x, y1, y2, 0, 1000);
end

end
