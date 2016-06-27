function [b] = UncertaintyAccuracy(field)
% plot RT as a function of |motionstrength|, to show that the pattern
% follows model predictions of uncertainty

global mypath;
nbins = 12;

% can try this also with all subjects
subjects = 1:27;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grandavg.rt     = nan(length(subjects), nbins);
grandavg.acc    = nan(length(subjects), nbins);
grandavg.b      = nan(length(subjects), 2);

for sj = subjects,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % bin
    [grandavg.rt(sj, :), grandavg.acc(sj, :)] = divideintobins(data.(field), data.correct, nbins);
    
    % project RT out of the pupil and vice versa
    switch field
        case 'rt'
            data.(field) = projectout(zscore(log(data.rt + 0.1)), zscore((data.decision_pupil)));
        case 'decision_pupil'
            data.(field) = projectout(zscore(data.decision_pupil), zscore(log(data.rt + 0.1)));
    end
    
    % also do logistic regression
    grandavg.b(sj, :) = glmfit(zscore(data.(field)), data.correct, ...
        'binomial','link','logit');
end

h = boundedline(1:nbins, 100* nanmean(grandavg.acc), ...
    100 * nanstd(grandavg.acc) ./ sqrt(length(subjects)), 'cmap', [0 0 0]);

% stats
xlim([0 nbins]); set(gca, 'xtick', [1 nbins/2 nbins], 'xminortick', 'off');
switch field
    case 'rt'
        ylim([45 90]); set(gca, 'ytick', [50:20:100]);
        xlabel('Reaction time');
        set(gca,  'xticklabel', {'fast', 'med', 'slow'});
    case 'decision_pupil'
        xlabel('Pupil response');
        ylim([68 80]); set(gca, 'ytick', [70:10:100]);
        set(gca,  'xticklabel', {'low', 'medium', 'high'});
end
ylabel('Accuracy (%)');
axis square;

b = grandavg.b;
end