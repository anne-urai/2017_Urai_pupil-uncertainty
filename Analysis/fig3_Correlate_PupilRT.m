% check to what extend the single-trial values of RT and pupil responses
% are correlated 

subjects = 1:27; clc
for sj = unique(subjects),
    
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data3_sj%02d.csv', sj));
    
    subplot(5,6,sj); plot(data.decision_pupil, data.rt, '.'); axis tight; box off;
    [grandavg.spearman(sj), grandavg.spearmanpval(sj)]  = corr(data.decision_pupil, data.rt, 'type', 'spearman');
    [grandavg.pearson(sj), grandavg.pearsonpval(sj)] = corr(data.decision_pupil, data.rt, 'type', 'pearson');
  %  ylim([0.5 3]);
   % disp(mean(data.decision_pupil));
end

% check across the group
suplabel('Pupil response', 'x'); suplabel('RT', 'y');

disp(mean(grandavg.pearson))
disp(min(grandavg.pearson))
disp(max(grandavg.pearson))

disp(length(find(grandavg.pearsonpval < 0.05)))


