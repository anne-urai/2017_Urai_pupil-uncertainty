% replicates figure 3
% Anne Urai, 2015

subplot(4,4,1); fig3a_PupilTimecourse;
subplot(4,4,3); fig3b_PupilUncertaintyTimecourse;

subplot(4,4,5); fig3c_PupilUncertaintyCorrelation;
subplot(4,4,7); fig3de_Uncertainty_Accuracy;

print(gcf, '-dpdf', sprintf('%s/Figures/figure3.pdf', mypath));

%% not in figure, but compute the correlation between RT and pupil values for each SJ
clc;
subjects = 1:27; 
for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    subplot(5,6,sj); plot(data.decision_pupil, data.rt, '.'); axis tight; box off;
    [grandavg.pearson(sj), grandavg.pearsonpval(sj)] = corr(data.decision_pupil, data.rt, 'type', 'pearson');
end

% check across the group
suplabel('Pupil response', 'x'); suplabel('RT', 'y');
disp(mean(grandavg.pearson))
disp(min(grandavg.pearson))
disp(max(grandavg.pearson))
disp(length(find(grandavg.pearsonpval < 0.05)))
