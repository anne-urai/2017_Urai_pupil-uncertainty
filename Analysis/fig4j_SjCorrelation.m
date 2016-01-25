function [] = fig4j_SjCorrelation(whichmodulator)
% correct - error = stim + resp --stim +- resp = 2stim
% error - correct = -stim + resp -resp - stim = -2stim
global mypath;

if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end

load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator));
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));

scatter(dat.response(:, 1), dat.response_pupil(:, 1), 10, mycolmap, 'filled');

[rho, pval] = corr(dat.response(:, 1), dat.response_pupil(:, 1), 'type', 'pearson');
max(abs(dat.response(:, 1)))
max(abs(dat.response_pupil(:, 1)))
if pval < 0.05,
    lsline;
end
axis square; box on;

text(0.5, -0.15, sprintf('r = %.3f', rho));
text(0.5, -0.2, sprintf('p < %.3f', pval));
ylim([-.25 .25]); set(gca, 'ytick', [-.2 0 0.2]);
xlim([-0.45 0.45]); set(gca, 'xtick', [-.4 0 0.4]);

bf10 = corrbf(rho,27);
text(0.5, -0.25, sprintf('bf10 = %.3f', bf10));

switch whichmodulator
    case 'pupil-rt'
        ylabel('Pupil');
    case 'rt-pupil'
        ylabel('RT');
        xlabel('Choice weight');
        
        % between subject stuff
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil-rt'));
        x = dat.response(:, 1);
        y1 = dat.response_pupil(:, 1);
        load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'rt'));
        y2 = dat.response_pupil(:, 1);
        
        [pval, deltaR] = permtest_correlation(x, y1, y2, 0, 10000);
        disp(pval);
end
set(gca, 'xcolor', 'k', 'ycolor', 'k');

end
