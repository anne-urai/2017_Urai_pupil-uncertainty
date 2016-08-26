function [] = PupilDecayFunction()
% ========================================================= %
% neuromodulatory decay
% ========================================================= %
global mypath;
nlags = 7; % just as in Fruend
hold on;

subjects = 1:27;
for sj = unique(subjects),
    data    = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % first, get the correlation at lag 0 between decision pupil and the
    % end of the trial
    neuromodDecay.sametrial(sj) = corr(data.decision_pupil, data.trialend_pupil, 'rows', 'complete', 'type', 'Spearman');    
    trls = 1:length(data.decision_pupil);
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

        % these look more or less the same
        neuromodDecay.pear(sj, lag) = corr(data.decision_pupil(trls), ...
            data.baseline_pupil(laggedtrls), 'type', 'Pearson', 'rows', 'complete');
    end
end

% split subjects based on their plain history weights

subplot(221);
hold on;
% plot([1 lag], [0 0], 'k', 'linewidth', 0.5);
plot(1:lag, neuromodDecay.pear', 'color', [0.7 0.7 0.7], 'linewidth', 0.5);
plot(1:lag, mean(neuromodDecay.pear), 'color', [0 0 0], 'linewidth', 1);
axisNotSoTight; axis square;
xlim([0 nlags+0.5]);
xlabel('Lags'); ylabel({'Spearman''s rho'});
set(gca, 'xtick', 1:nlags);
% ylim([-0.08 0.12]);
set(gca, 'xcolor', 'k', 'ycolor', 'k');

