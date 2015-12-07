
% ========================================================= %
% panel E: neuromodulatory decay
% ========================================================= %
clf; clear; clc;
nlags = 5;

subjects = 1:27;
for sj = unique(subjects),
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    trls = 1:size(data, 1);
    
    % first, get the correlation at lag 0 between decision pupil and the
    % end of the trial
    tmp = corrcoef(data.decision_pupil', data.trialend_pupil', 'rows', 'complete');
    neuromodDecay.sametrial(sj) = tmp(2);

    for lag = 1:nlags,
        
        % with this selection, take the trials after that
        laggedtrls = trls+lag;
        
        % exclude trials at the end of the block
        if any(laggedtrls > size(data, 1)),
            trls(laggedtrls > size(data, 1)) = [];
            laggedtrls(laggedtrls > size(data, 1)) = [];
        end
        % remove trials that dont match in block nr
        noblockmatch = data.blocknr(laggedtrls) ~= data.blocknr(trls);
        laggedtrls(noblockmatch) = [];
        trls(noblockmatch) = [];
        
        % subtract correct the baseline on tLagged by the baseline on t0
        % baselineCorr = projectout(data.baseline_pupil(laggedtrls), data.baseline_pupil(trls));
        % subtract or project out?
        baselineCorr = data.baseline_pupil(laggedtrls) - data.baseline_pupil(trls);
        
        % these look more or less the same
        %  neuromodDecay.spear(sj, lag) = corr(baselineCorr, data.decision_pupil(trls), 'type', 'Pearson');
        tmp = corrcoef(baselineCorr, data.decision_pupil(trls));
        neuromodDecay.pear(sj, lag) = tmp(2);
        
    end
end

subplot(4,4,1); hold on;
errorbar(0, nanmean(neuromodDecay.sametrial), (nanstd(neuromodDecay.sametrial) ./ sqrt(length(subjects))), ...
    '.k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'markersize', 8);
b1 = boundedline(1:lag, mean(neuromodDecay.pear(:, :, 1)), std(neuromodDecay.pear(:, :, 1)) ./ sqrt(length(subjects)), 'cmap', [0.5 0.5 0.5]);
%b2 = boundedline(1:lag, mean(neuromodDecay.spear), std(neuromodDecay.spear) ./ sqrt(length(subjects)), 'cmap', [0.3 0.3 0.3]);
axis tight;  xlim([-1 lag+0.5]);
xlabel('Lags'); ylabel({'Pearson''s rho'});
set(gca, 'xtick', 0:nlags);
ylim([0.4 0.7]);

%lh = legend([b1 b2], {'Pearson', 'Spearman'}, 'Location', 'NorthEast'); legend boxoff;
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4b_neuroModDecay2.pdf'));
