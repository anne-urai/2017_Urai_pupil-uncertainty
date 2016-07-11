% figure 4 overview
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
subplot(4,4,1); plotBetasSwarm([dat.response_pupil(:, 1) ...
    dat.stimulus_pupil(:, 1)  dat.response_rt(:, 1) dat.stimulus_rt(:, 1)], ...
    [0 0 0; 0 0 0; 0 0 0; 0 0 0]);
set(gca, 'xtick', 1:4, 'xticklabel', ...
    {'Pupil * choice', 'Pupil * stimulus', 'RT * choice', 'RT * stimulus'}, ...
    'xticklabelrotation', -30);

subplot(5,5,3); SjCorrelation('pupil', 'response');
subplot(5,5,4); SjCorrelation('rt', 'response');

% show median split for correlation stuff
subplot(4,6, 13); MedianSplit('pupil', 'response'); 
ylim([-0.4 0.15]); ylabel('Pupil * choice');
subplot(4,6,14); MedianSplit('rt', 'response'); 
ylim([-0.4 0.15]); set(gca, 'yaxislocation', 'right');
ylabel('RT * choice');

print(gcf, '-dpdf', sprintf('%s/Figures/figure6.pdf', mypath));

