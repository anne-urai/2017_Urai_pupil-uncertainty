% does decision pupil have an effect on the probability to stay vs switch responses on the next trial?
% Anne Urai, UKE, 2015

%% show switching pattern per pupil bin
clc; clear; close all;

subjects = 1:27;
nbins = 20; % for display purposes only, low-medium-high

for sj = unique(subjects),
    data = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data_sj%02d.csv', sj));
    
    % predict response
    resp = data.resp; resp(resp == -1) = 0;
    
    % don't use trials that are at the end of a block (since reasonably,
    % they won't predict the first response after the break/end of the day)
    trlDif = [diff(data.trialnr); NaN];
    resp(trlDif < 1) = nan;
    resp(trlDif > 1) = nan;
    resp(find(trlDif > 1) + 1) = nan;
    
    
    % run with only 1 lag
    designM = [zscore(data.motionstrength) ...
        circshift(data.resp, 1) ...
        circshift(zscore(data.decision_pupil), 1) ...
        circshift(data.resp, 1) .* circshift(zscore(data.decision_pupil), 1) ...
        ];
    
    xlabs = {'bias', 'stim', 'resp-1', 'pup-1', 'respPup-1'};
    
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
        grandavg.betas(sj, c, :) = mdl.Coefficients.Estimate;
        grandavg.pval(sj, c, :)  = mdl.Coefficients.pValue;
    end
end

%% show the result of the full logistic regression model
clf;
colors = linspecer(3); %colors = colors([2 3 1], :);

for sp = 1:2,
    
    if sp == 1,
        what2plot = 1:2; % bias stimswitch
    elseif sp == 2,
        what2plot = 3:5; % resp pupil interaction
    end
    
    subplot(4,4,sp); hold on;
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
    set(gca, 'xtick', 1:length(what2plot), 'xticklabel', []);
    
    if sp == 1,
        ylabel('Beta weights (a.u.)');
    elseif sp == 2,
        set(gca,'yaxislocation','right');
    end
    box off;
end

% now, do the same but make sure we have red and bar graphs
colors = linspecer(3); colors = colors([2 3 1], :);

for sp = 1:2,
    
    if sp == 1,
        what2plot = 1:2; % bias stimswitch
    elseif sp == 2,
        what2plot = 3:5; % resp pupil interaction
    end
    
    subplot(4,4,sp+4); hold on;
    
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
    % h = sigstar(arrayfun(@repmat, x, ones(1, n), 2*ones(1, n), 'UniformOutput', 0), pval);
    
    xlims = get(gca, 'xlim'); set(gca, 'xlim', [min(xlims)-0.5 max(xlim)+0.5]);
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    set(gca, 'xtick', 1.5:2:length(what2plot)*2, 'xticklabel', xlabs(what2plot));
    
    if sp == 1,
        ylabel('Beta weights (a.u.)');
    elseif sp == 2,
        set(gca,'yaxislocation','right');
    end
    box off;
end


print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4_sequentialDependencies.pdf'));

