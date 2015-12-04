%% see if the extra parameters are warranted by the AIC

aic = @(logL, numParam) -2*(logL)+2*(numParam);

subjects = 1:20;
for sj = subjects,
    
    % first get the model without pupil
    load(sprintf('~/Data/UvA_pupil/sim_backup/2ifc_sj%02d.txtresults.mat', sj));
    
    dat.aic(sj, 1) = aic(model_nohist.loglikelihood, length(model_nohist.w0));
    dat.aic(sj, 2) = aic(model_w_hist.loglikelihood, length(model_w_hist.w0));
    
    % then include the pupil term
    load(sprintf('~/Data/UvA_pupil/sim_backup/2ifc_pupil_sj%02d.txtresults.mat', sj));
    dat.aic(sj, 3) = aic(model_w_hist.loglikelihood, length(model_w_hist.w0));
end

% normalize by the AIC in the no history model
dat.aicNorm = bsxfun(@minus, dat.aic, dat.aic(:, 1));

subplot(221);
boxplot(dat.aicNorm(:, 2:end));

set(gca, 'xtick', 1:2, 'xticklabel', {'history model', 'history + pupil'});
ylabel('AIC, normalized by no history model');

% test whether the pupil model has a lower aic
[h, p, ci, stats] = ttest(dat.aicNorm(:, 2), dat.aicNorm(:, 3), 'tail', 'right')