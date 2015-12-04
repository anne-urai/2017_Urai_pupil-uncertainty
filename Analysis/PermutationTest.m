clear; clf;

subjects = 1:20;
for sj = subjects,
    
    % first get the model without pupil
    load(sprintf('~/Data/UvA_pupil/sim_backup/2ifc_sj%02d.txtresults.mat', sj));
    
    subplot(5,6,sj);
    histogram(permutation_wh(:, 1), 'EdgeColor', 'none'); hold on;
    axis tight; box off; 
    plot([model_w_hist.loglikelihood model_w_hist.loglikelihood], [0 max(get(gca, 'ylim'))])
end

suplabel('Loglikelihood of model with history coupling');


%% do the same thing but with pupil
clear; clf;

subjects = 1:20;
for sj = subjects,
    
    % first get the model without pupil
    load(sprintf('~/Data/UvA_pupil/sim_backup/2ifc_sj%02d.txtresults.mat', sj));
    
    subplot(5,6,sj);
    histogram(permutation_wh(:, 1), 'EdgeColor', 'none'); hold on;
    axis tight; box off; 
    plot([model_w_hist.loglikelihood model_w_hist.loglikelihood], [0 max(get(gca, 'ylim'))])
    
end

suplabel('Loglikelihood of model with history coupling');

