
% reproduces 
global mypath;

close;
figure;

subplot(442); historyContribution;

% take P10 as an example
subplot(443); FruendKernels('plain', 'response');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

subplot(444); decisionStrategies('plain');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
subplot(4,4,9); plotBetas([dat.response_pupil(:, 1) ...
    dat.stimulus_pupil(:, 1)  dat.response_rt(:, 1) dat.stimulus_rt(:, 1)], ...
    0.5*ones(4,3));

% paired stats
[~, pval, ~, stat] = ttest(dat.response_pupil(:, 1), dat.stimulus_pupil(:, 1));
mysigstar(gca, [1 2], -0.08, pval);
[~, pval, ~, stat] = ttest(dat.response_rt(:, 1), dat.stimulus_rt(:, 1));
mysigstar(gca, [3 4], -0.08, pval);

xlim([0.5 4.5]);
set(gca, 'xtick', 1:4, 'xticklabel', ...
    {'Pupil x choice', 'Pupil x stimulus', 'RT x choice', 'RT x stimulus'}, ...
    'xticklabelrotation', -30);

subplot(4,4,10); SjCorrelation('pupil', 'response');
subplot(4,4,11); SjCorrelation('rt', 'response');

% show median split for correlation stuff
subplot(4,9,35); MedianSplit('pupil', 'response'); 
ylim([-0.4 0.15]); ylabel('Pupil x choice weight');
subplot(4,9,36); MedianSplit('rt', 'response'); 
ylim([-0.4 0.15]); set(gca, 'yaxislocation', 'right');
ylabel('RTx choice weight');

print(gcf, '-dpdf', sprintf('%s/Figures/Figure5.pdf', mypath));


