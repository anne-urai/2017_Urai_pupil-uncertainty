% rather than waiting for the history model to finish in Python, run it
% myself

%% show switching pattern per pupil bin
clc; clear; close all;

subjects = 1:27;
grandavg.betas = nan(27, 3, 7);
grandavg.pval = nan(27, 3, 7);

for sj = unique(subjects),
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    % get the variables we need
    resp = data.resp; resp(resp == -1) = 0; % predict response identity
    motionstrength = zscore(data.motionstrength);
    prevPupil      = circshift(zscore(data.decision_pupil), 1);
    prevResp       = circshift(data.resp, 1);
    prevStim       = circshift(data.stim, 1);
    
    % make design matrix - intercept will be added automatically
    designM = [motionstrength prevResp prevResp.*prevPupil];
    xlabs = {'bias', 'stim', 'prevResp','prevResp*prevPupil'};
    
    % don't use trials that are at the beginning of each block
    trlDif = [0; diff(data.trialnr)];
    
    removeTrls = false(size(trlDif));
    removeTrls(trlDif < 1) = true;
    removeTrls(trlDif > 1) = true;
    removeTrls(find(trlDif > 1) + 1) = true;
    
    % in these trials, the history wont be able to predict the response
    designM(removeTrls==1, 2:end) = 0;
     
    % do this separately for error, correct and all trials together
    cors = [0 1];
    for c = 1:3,
        correctPrev = circshift(data.correct, 1);
        if c == 3,
            trls = 1:length(resp);
        else
            trls = find(correctPrev == cors(c));
        end

        mdl  = fitglm(designM(trls, :), resp(trls), ...
            'distr', 'binomial', 'link', 'logit');
        grandavg.betas(sj, c, 1:length(mdl.Coefficients.Estimate)) = mdl.Coefficients.Estimate;
        grandavg.pval(sj, c, 1:length(mdl.Coefficients.pValue))  = mdl.Coefficients.pValue;
        
        grandavg.dev(sj, 3) = mdl.Deviance;
        
        if c == 3, % only on all trials
            
            % get the deviance for the null model without history
            mdl  = fitglm(designM(trls, 1), resp(trls), ...
                'distr', 'binomial', 'link', 'logit');
            grandavg.dev(sj, 1) = mdl.Deviance;
            
            % model with  history
            mdl  = fitglm(designM(trls, 1:2), resp(trls), ...
                'distr', 'binomial', 'link', 'logit');
            grandavg.dev(sj, 2) = mdl.Deviance;
            
        end
    end
end

savefast('~/Data/pupilUncertainty/GrandAverage/historyMatlab.mat', 'grandavg', 'xlabs');

%% deviance statistics between the different models

dev.chisqHist = grandavg.dev(:, 1) - grandavg.dev(:, 2);
dev.histPval = 1 - chi2cdf(dev.chisqHist,1);

dev.chisqPup = grandavg.dev(:, 2) - grandavg.dev(:, 3);
dev.pupilPval = 1 - chi2cdf(dev.chisqPup,1);

%% show the result of the full logistic regression model
clf;
colors = linspecer(3); % colors = colors([2 3 1], :);

for sp = 1:2,
    
    if sp == 1,
        what2plot = 1:2; % bias stimswitch
    elseif sp == 2,
        what2plot = 3:4; % resp, stim,  pupil interaction
    end
    
    subplot(2,2,sp); hold on;
    % start with a graph for all the trials combined
    bar(1:length(what2plot), squeeze(mean(grandavg.betas(:, 3, what2plot))), ...
        'facecolor', colors(1, :), 'edgecolor', 'none', 'barwidth', 0.5);
    h = ploterr(1:length(what2plot), squeeze(mean(grandavg.betas(:, 3, what2plot))), ...
        [], squeeze(std(grandavg.betas(:, 3, what2plot))) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    axis tight; axis square;
    
    clear pval
    for n = 1:length(what2plot),
        [~, pval(n)] = permtest(squeeze(grandavg.betas(:, 3, what2plot(n))));
    end
    h = sigstar(arrayfun(@repmat, 1:n, ones(1, n), 2*ones(1, n), 'UniformOutput', 0), pval);
    
    xlims = get(gca, 'xlim'); set(gca, 'xlim', [min(xlims)-0.5 max(xlim)+0.5]);
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    set(gca, 'xtick', 1:length(what2plot), 'xticklabel', xlabs(what2plot), 'xticklabelrotation', -30);
    
    if sp == 1,
        ylabel('Beta weights (a.u.)');
    elseif sp == 2,
        set(gca,'yaxislocation','right');
    end
    box off;
end

% now, do the same but separately for correct and error
colors = linspecer(3); colors = colors([2 3 1], :);

for sp = 1:2,
    
    if sp == 1,
        what2plot = 1:2; % bias stimswitch
    elseif sp == 2,
        what2plot = 3:4; % resp pupil interaction
    end
    
    subplot(2,2,sp+2); hold on;
    
    % correct and error
    for c = 1:2,
        
        switch c
            case 1
                x = 1:2:length(what2plot)*2;
            case 2;
                x = 2:2:length(what2plot)*2;
        end
        
        bar(x, squeeze(mean(grandavg.betas(:, c, what2plot))), ...
            'facecolor', colors(c, :), 'edgecolor', 'none', 'barwidth', 0.4);
        h = ploterr(x, squeeze(mean(grandavg.betas(:, c, what2plot))), ...
            [], squeeze(std(grandavg.betas(:, c, what2plot))) / sqrt(length(subjects)), ...
            'k.', 'hhxy', 0.0000000000000001);
        set(h(1), 'marker', 'none');
        axis tight; axis square;
        
        % within error and correct, compute significance across the group
        clear pval
        for n = 1:length(what2plot),
            [~, pval(n)] = permtest(squeeze(grandavg.betas(:, c, what2plot(n))));
        end
        h = sigstar(arrayfun(@repmat, x, ones(1, n), 2*ones(1, n), 'UniformOutput', 0), pval);
    end
    
    % also add sigstars to compare each regressor between error and correct
    clear pval
    for n = 1:length(what2plot),
        [~, pval(n)] = permtest(squeeze(grandavg.betas(:, 1, what2plot(n))), squeeze(grandavg.betas(:, 2, what2plot(n))));
        sigstar({[x(n)-1, x(n)]}, pval(n));
    end
    
    xlims = get(gca, 'xlim'); set(gca, 'xlim', [min(xlims)-0.5 max(xlim)+0.5]);
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    set(gca, 'xtick', 1.5:2:length(what2plot)*2, 'xticklabel', xlabs(what2plot), 'xticklabelrotation', -30);
    
    if sp == 1,
        ylabel('Beta weights (a.u.)');
    elseif sp == 2,
        set(gca,'yaxislocation','right');
    end
    box off;
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4_sequentialDependencies.pdf'));

%% scatter plot to see if pupil weight depends on history weights
figure;
colors = linspecer(3); 
for i = 1:3,

    % original weight and interaction term
    mainIdx = find(strcmp(xlabs, 'prevResp')==1);
    intIdx = find(strcmp(xlabs, 'prevResp*prevPupil')==1);
    switch i
        case 1
            mainW = grandavg.betas(:, 3, mainIdx);
            intW  = grandavg.betas(:, 3, intIdx);
        case 2
            mainW = grandavg.betas(:, 1, mainIdx);
            intW = grandavg.betas(:, 1, intIdx);
        case 3
            mainW = grandavg.betas(:, 2, mainIdx);
            intW = grandavg.betas(:, 2, intIdx);
    end
    
    subplot(4,4,i+12);
    plot(mainW, intW, 'o', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'w');
    
    [rho, pval] = corr(mainW, intW, 'type', 'spearman');
    if pval < 0.05,
        ls = lsline; set(ls, 'color', [0.7 0.7 0.7]);
    end
    disp([rho pval]);
    box off;

end


