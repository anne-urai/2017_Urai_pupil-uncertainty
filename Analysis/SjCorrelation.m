function [] = SjCorrelation(whichmodulator, whichweight, whichFile)
% plot correlation between subjects

if ~exist('whichweight', 'var'); whichweight = 'response'; end
if ~exist('whichFile', 'var'); whichFile = 'pupil+rt'; end

global mypath;
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichFile));
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));

hold on;
% axes
plot([-1 1], [0 0], 'color', [0.5 0.5 0.5], 'linewidth', 0.2);
plot([0 0 ], [-1 1], 'color', [0.5 0.5 0.5], 'linewidth', 0.2);

% or without errorbars
scatter(dat.response(:, 1), dat.([whichweight '_' whichmodulator])(:, 1), 10, mycolmap, 'filled');

[rho, pval] = corr(dat.response(:, 1), dat.([whichweight '_' whichmodulator])(:, 1), 'type', 'pearson');
%if pval < 0.05,
lsline;
%end

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

title(sprintf('r = %.3f, p = %.3f', rho, pval));
ylim([-0.45 0.25]); set(gca, 'ytick', [-.4:0.2:0.4]);
xlim([-0.5 0.5]); set(gca, 'xtick', [-.4 0 0.4]);

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
        
        % importantly, in this case they return the same value
end
axis square; box on;

end
