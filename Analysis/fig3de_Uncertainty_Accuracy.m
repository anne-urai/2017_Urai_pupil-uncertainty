function [grandavg] = fig3de_Uncertainty_Accuracy(nbins)
% plots uncertainty by accuracy both for the modelfits and the pupil

if ~exist('nbins', 'var'), nbins = 100; end
close all;

subjects      = 1:27;
fitIndividual = false;

for sj = unique(subjects),
    
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    % divide into bins
    [ grandavg.pup(sj, :), grandavg.acc(sj, :), stdx, stdy] = ...
        divideintobins(data.decision_pupil, data.correct, nbins);
    
    % do a logistic regression
    x = [zscore(abs(data.motionstrength)) zscore(data.rt) zscore(data.decision_pupil) ];
    
    x(data.trialnr < 2, :) = nan;
    
    [b,dev,stats] = glmfit(x, data.correct, ...
        'binomial','link','logit');
    % save for later
    grandavg.b(sj, :) = b;
    
    if fitIndividual,
        subplot(5,6,sj);
        yfit = glmval(b, x,'logit');
        hold on;
        plot(grandavg.pup(sj, :), grandavg.acc(sj,:), 'b.');
        axis tight; ylim([0.5 1]); xlims = get(gca, 'xlim');
        plot(x,yfit,'k');
        xlim(xlims);
        title(sprintf('P%02d', sj));
    end
    
end

% now, plot
cols = linspecer(3);
figure;
subplot(441);

grandavg.acc = 100 * grandavg.acc;
boundedline(1:nbins, nanmean(grandavg.acc), ...
    nanstd(grandavg.acc) / sqrt(length(subjects)), 'cmap', [0 0 0]);

xlabel('Pupil response');
ylabel({'Percent correct'});
axis tight; axis square;
set(gca, 'ytick', [65 75 85]);
set(gca, 'xtick', [1  nbins/2 nbins], 'xticklabel', {'low', 'medium', 'high'});
% offsetAxes;
box off; ylim([65 85]);

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig3d_pupil_accuracy.pdf'));

%% fit logistic slope on high vs low pupil bins
subjects = 1:27;
grandavg.betas = nan(length(subjects), 2, 2);

for sj = unique(subjects),
    
    disp(sj);
    % get all the data
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    % divide into low and high pupil bins
    puptrls{1} = find(data.decision_pupil < quantile(data.decision_pupil, 0.5));
    puptrls{2} = find(data.decision_pupil > quantile(data.decision_pupil, 0.5));
    
    % clear puptrls; % puptrls{1} = 1:height(data);
    
    for b = 1:length(puptrls),
        thisdat = data(puptrls{b}, :);
        
        % fit on all the trials, dont average beforehand!
        x = thisdat.motionstrength;
        y = thisdat.resp; y(y==-1) = 0;
        
        % sort both
        [x, idx] = sort(x);
        y = y(idx);
        
        [beta,dev,stats] = glmfit(x, [y ones(size(y))], ...
            'binomial','link','logit');
        newx = -3:0.01:3;
        yfit = glmval(beta, newx,'logit');
        % plot(newx, yfit,'-', 'color', cols(b, :));
        
        % also save the curve itself
        grandavg.betas(sj, b, :) = beta;
        grandavg.yfit(sj, b, :) = yfit;
        
        % also get the individual d prime
        % use 0 (for weaker motion answer) instead of -1
        resp = thisdat.resp; resp(resp==-1) = 0;
        stim = thisdat.stim; stim(stim==-1) = 0;
        
        hit  = length(find(resp(find(stim == 1))==1)) / length(resp(find(stim == 1)));
        fa   = length(find(resp(find(stim == 0))==1)) / length(resp(find(stim == 0)));
        
        grandavg.dprime(sj, b) = norminv(hit) - norminv(fa);
        
        % see if there is a difference in the actual stimulus as well
        grandavg.difficulty(sj, b) = mean(abs(thisdat.motionstrength));
    end
end

% make sure I actually plot the bar graph here to reproduce figure 3!
subplot(4,5,3);
h = ploterr(1:2, squeeze(nanmean(grandavg.betas(:, :, 2))), [], ...
    squeeze(std(grandavg.betas(:, :, 2))) / sqrt(length(subjects)), 'k.', 'hhxy', 0.00001);
set(h(1), 'marker', 'none');

hold on;
bar(1, mean(grandavg.betas(:, 1, 2)), 'facecolor', [0.6 0.6 0.6], 'edgecolor', 'none', 'barwidth', 0.4);
bar(2, mean(grandavg.betas(:, 2, 2)), 'facecolor', [0.4 0.4 0.4], 'edgecolor', 'none', 'barwidth', 0.4);
% do ttest on regression coefficients
[~, pval(3), ~, stat] = ttest(grandavg.betas(:, 1, 1), grandavg.betas(:, 1, 2));
ylim([0.9 1.5]); xlim([0.5 2.5]); box off;
s1 = sigstar({[1 2]}, pval(3), 0); ylim([0.9 1.5]); set(gca, 'ytick', [1 1.5]);
xlabel('Pupil response'); set(gca, 'xtick', 1:2, 'xticklabel', {'low', 'high'});
ylabel('Perceptual sensitivity');

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig3e_pupil_sensitivity.pdf'));
savefast('~/Data/pupilUncertainty/GrandAverage/grandavg_logistic_bypupil.mat', 'grandavg', 'subjects', 'cols');

%%
% !!! rather than doing a permutation test on the logistic slope
% coefficients, run a proper mixed effects logistic regression
% fixed effects: slope
if 0,
    clear; clc;
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_allsj.csv'));
    
    % normalize the pupil response for each subject
    % this way, we can interpret the coefficient of the pupil fixed effect in
    % units of standard deviation
    for sj = unique(data.subjnr)',
        data.decision_pupil(sj==data.subjnr) = zscore(data.decision_pupil(sj==data.subjnr));
    end
    
    data.absmotion = abs(data.motionstrength);
    
    mdl = fitglme(data, 'correct ~ 1 + decision_pupil + absmotion + (1|subjnr) + (1|sessionnr) ', ...
        'Distribution', 'Binomial', 'Link', 'Logit');
    disp(mdl);
end

%% compute ROC AUC based on pupil
clear all; clc; close all;

subjects = 1:27;
for sj = unique(subjects),
    
    % get all the data
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    out = rocAnalysis(data.decision_pupil(data.correct==1), ...
        data.decision_pupil(data.correct==0), 0, 1);
    
    grandavg.roc(sj)    = out.i;
    grandavg.pval(sj)   = out.p;
end
savefast('~/Data/pupilUncertainty/GrandAverage/pupilCorrectnessROC.mat', 'grandavg');

load('~/Data/pupilUncertainty/GrandAverage/pupilCorrectnessROC.mat');
% test across the group, roc auc values are normally distributed
[h, pval ci, stats] = ttest(grandavg.roc, 0.5, 'tail', 'both')

end
