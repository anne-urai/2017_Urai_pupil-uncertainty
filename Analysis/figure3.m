% replicates figure 3
% Anne Urai, 2015

close all; 
global mypath;

subplot(4,4,1); PupilTimecourse(0);
subplot(4,4,3); PupilUncertaintyTimecourse;

% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 9);
% error vs correct
subplot(445); b = Uncertainty_byErrorCorrect('decision_pupil');
cla; plotBetasSwarm(b, colors([1 2], :));
set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'});

% other metrics of uncertainty
subplot(4,4,6); b = UncertaintyAccuracy('decision_pupil');
cla; plotBetasSwarm(b(:, 2), [0 0 0]);
set(gca, 'xtick', 1, 'xticklabel', []);

% psychometric functions
subplot(4,4,7); b = PsychFuncs_byUncertainty('decision_pupil');
cla; plotBetasSwarm(1./b, [0.7 0.7 0.7; 0.2 0.2 0.2]);
set(gca, 'xtick', [1 2], 'xticklabel', {'low', 'high'});
xlabel('Pupil response'); ylabel('Sensitivity (a.u.)');
ylim([0 1]);

% ensure same axes proportions
ax = findobj(gcf, 'type', 'axes');
for a = 1:length(ax),
    pos = get(ax(a), 'position'); pos(4) = 0.15; set(ax(a), 'position', pos);
end
print(gcf, '-dpdf', sprintf('%s/Figures/Figure3.pdf', mypath));
% note: several of these panels appear in Figure S3 in the paper.

%% compute the correlation between RT and pupil values for each SJ
clc;
subjects = 1:27; 
for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data.rt = zscore(log(data.rt + 0.1));
    
    %subplot(5,6,sj); plot(data.decision_pupil, data.rt, '.'); axis tight; box off;
    [grandavg.pearson(sj), grandavg.pearsonpval(sj)] = corr(data.decision_pupil, data.rtNorm, 'type', 'pearson');
end
% check across the group
fprintf('average r: %.3f range: %.3f to %.3f, significant in %d out of %d observers \n', ...
    mean(grandavg.pearson), min(grandavg.pearson), max(grandavg.pearson), length(find(grandavg.pearsonpval < 0.05)), 27)
