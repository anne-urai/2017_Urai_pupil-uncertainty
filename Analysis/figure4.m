% replicates figure 4
% Anne Urai, 2015

close all; figure;
global mypath;

subplot(4,4,1); PupilTimecourse(0);
subplot(4,4,3); PupilUncertaintyTimecourse;
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));

% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 9);

% error vs correct
subplot(445); b = Uncertainty_byErrorCorrect('decision_pupil');
subplot(4,4,6); scatterHistDiff(squeeze(b(:, 1)), squeeze(b(:, 2)), [], [], mycolmap);
xlabel('Correct trials beta weight (a.u.)'); ylabel('Error trials beta weight (a.u.)');

%subplot(4,7,10); plotBetasSwarm(b(:, :, 2), colors([1 2], :));
%set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'});

% other metrics of uncertainty
subplot(4,4,9); b = UncertaintyAccuracy('decision_pupil');
%subplot(4,4,10); histogram(b(:, 2));

% show betas
subplot(4,7,17); plotBetasSwarm(b(:, 2), [0 0 0]);
set(gca, 'xtick', 1, 'xticklabel', []);
%ylim([-0.6 0]);

% psychometric functions
subplot(4,4,13); b = PsychFuncs_byUncertainty('decision_pupil');
subplot(4,4,14); scatterHistDiff(squeeze(b(:, 1)), squeeze(b(:, 2)), ...
   [], [], mycolmap);

%subplot(4,7,24); plotBetasSwarm(b(:, :, 2), [0.7 0.7 0.7; 0.2 0.2 0.2]);
xlabel('Low pupil threshold'); ylabel('High pupil threshold');
%set(gca, 'xtick', [1 2], 'xticklabel', {'low', 'high'});
%xlabel('Pupil response'); ylabel('Threshold (a.u.)');
%ylim([0 7]);

% ensure same axes proportions
ax = findobj(gcf, 'type', 'axes');
for a = 1:length(ax),
    pos = get(ax(a), 'position'); pos(4) = 0.15; set(ax(a), 'position', pos);
end
print(gcf, '-dpdf', sprintf('%s/Figures/figure4.pdf', mypath));

%% not in figure, but compute the correlation between RT and pupil values for each SJ
clc;
subjects = 1:27; 
for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data.rt = zscore(log(data.rt + 0.1));
    
    %subplot(5,6,sj); plot(data.decision_pupil, data.rt, '.'); axis tight; box off;
    [grandavg.pearson(sj), grandavg.pearsonpval(sj)] = corr(data.decision_pupil, data.rt, 'type', 'pearson');
end

% check across the group
fprintf('pupil RT correlation, mean %f, min %f, max %f. nr of subjects with significant correlation: %d. \n', ...
    mean(grandavg.pearson), min(grandavg.pearson), max(grandavg.pearson), length(find(grandavg.pearsonpval < 0.05)))

%% rsquare for regresson models
load(sprintf('%s/Data/GrandAverage/pupilRegressionBetas.mat', mypath));
fprintf('R squared across regression samples: mean %f, std %d \n', mean(grandavg.rsq(:)), std(grandavg.rsq(:)));
