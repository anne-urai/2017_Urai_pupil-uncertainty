mypath = '/Users/anne/Data/pupilUncertainty_FigShare';
grandavg = table;

for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data(data.sessionnr > 1, :);
    
    % get the variables we need
    data.motionstrength = zscore(data.motionstrength);
    data.prevPupil      = circshift(zscore(data.decision_pupil), 1);
    data.prevRT         = circshift(zscore(log(data.rt + 0.1)), 1);
    data.prevUnc        = circshift(zscore(data.uncertainty), 1);
    data.prevResp       = circshift(data.resp, 1);
    data.prevStim       = circshift(data.stim, 1);
    data.resp(data.resp == -1) = 0; % predict response identity
    
    % don't use trials that are at the beginning of each block
    trlDif = [0; diff(data.trialnr)];
    removeTrls = false(size(trlDif));
    removeTrls(trlDif < 1) = true;
    removeTrls(trlDif > 1) = true;
    removeTrls(find(trlDif > 1) + 1) = true;
    data.resp(removeTrls) = NaN;
    
    % fit model
    modelSpec = ['resp ~ 1 + motionstrength + prevResp + prevStim' ...
        ' + prevPupil + prevPupil*prevResp + prevPupil*prevStim' ...
        ' + prevRT + prevRT*prevResp + prevRT*prevStim' ...
        ' + prevUnc + prevUnc*prevResp + prevUnc*prevStim'];
    mdl  = fitglm(data, modelSpec, 'distr', 'binomial', 'link', 'logit');
    
    % reformat
    coefnames = mdl.CoefficientNames;
    coefnames{1} = 'intercept';
    coefnames = regexprep(coefnames, ':', 'X');
    
    if sj == 1,
        grandavg = array2table(mdl.Coefficients.Estimate', 'VariableNames', coefnames);
    else
        grandavg = [grandavg; array2table(mdl.Coefficients.Estimate', 'VariableNames', coefnames)];
    end
end

%% show this
grandavg(:, 2) = []; % remove motionstrength, strongest beta

subplot(3,3,[1 3]);
plotBetasSwarm(table2array(grandavg));
set(gca, 'xtick', 1:width(grandavg), ...
    'xticklabel', grandavg.Properties.VariableNames, 'xticklabelrotation', -30);

% correlate all the 3 uncertainty * choice regressors with choice regressor
corrvars = {'prevPupil_prevResp', 'prevRT_prevResp', 'prevUnc_prevResp'};
for c = 1:length(corrvars),
    subplot(3,3,c+6);
    scatter(grandavg.prevResp, grandavg.(corrvars{c}), 'filled');
    xlim([-.5 0.5]); ylim([-.4 .8]);
    [rho, pval] = corr(grandavg.prevResp, grandavg.(corrvars{c}), 'type', 'spearman');
    if pval < 0.05,
        lsline;
    end
    ylabel(corrvars{c});
end
suplabel('prevResp', 'x');

print(gcf, '-dpdf', sprintf('%s/Figures/SwitchingThreeRegressors.pdf', mypath));

%%
figure;
grandavgResp = grandavg(:, [5 7 8 9]); % remove intercept, prevPupil, prevRT, prevUnc
corrplot(grandavgResp, grandavgResp.Properties.VariableNames);
set(findall(gcf,'type','text'),'FontSize',10);
print(gcf, '-dpdf', sprintf('%s/Figures/SwitchingThreeRegressorsMassiveCorrplotResp.pdf', mypath));

figure;
grandavgStim = grandavg(:, [6 10 11 12]); % remove intercept, prevPupil, prevRT, prevUnc
corrplot(grandavgStim, grandavgStim.Properties.VariableNames);
set(findall(gcf,'type','text'),'FontSize',10);
print(gcf, '-dpdf', sprintf('%s/Figures/SwitchingThreeRegressorsMassiveCorrplotStim.pdf', mypath));
