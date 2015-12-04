%% 1. dependence on pupil-based uncertainty, all trials
% does this pattern increase with higher uncertainty on previous trial?
clear; clc; close;
nbins = 3;
subjects = 1:27;

whichMod = 'decision_pupil';
whichLogisticParam = 2; % 1 = bias, 2 = slope
spcnt = 0;
regressoutRT = false;
for whichLag = 1:7,
    clearvars -except whichMod spcnt sj subjects nbins whichLag grandavg whichLogisticParam regressoutRT;
    
    for sj = unique(subjects),
        
        data = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data2_sj%02d.csv', sj));
        switch whichMod
            case 'rt'
                data.rt = log(data.rt);
            case 'decision_pupil'
                if regressoutRT,
                    data.decision_pupil = projectout(data.decision_pupil, data.rt);
                end
        end
        
        % outcome vector need to be 0 1 for logistic regression
        data.resp(data.resp == -1) = 0;
        
        % remove the trials at the end of each block
        endOfBlockTrls = find(data.trialnr == 50);
        
        % get an overall logistic fit
        [b, dev, stats] = glmfit(data.motionstrength, data.resp, ...
            'binomial','link','logit');
        grandavg.overallLogistic(sj, :) = b;
        
        % previous response
        resps = [0 1];
        for r = 1:2,
            
            if nbins > 2,
                uncQs = quantile(data.(whichMod)(data.resp == resps(r)), nbins - 1);
            elseif nbins == 2,
                uncQs = median(data.(whichMod)(data.resp == resps(r)));
            end
            
            % uncertainty bins
            for u = 1:nbins,
                
                switch u
                    case 1
                        trls = find(data.resp == resps(r) & data.(whichMod) < uncQs(u));
                    case nbins
                        trls = find(data.resp == resps(r) & data.(whichMod) > uncQs(u-1));
                    otherwise
                        trls = find(data.resp == resps(r) & ...
                            data.(whichMod) > uncQs(u-1) & data.(whichMod) < uncQs(u));
                end
                
                % only use correct trials!
                % trls = intersect(trls, find(data.correct == 1));
                
                % with this selection, take the trials after that
                usetrls = trls+whichLag;
                % exclude trials at the end of the block
                if any(usetrls > size(data, 1)),
                    trls(usetrls > size(data, 1)) = [];
                    usetrls(usetrls > size(data, 1)) = [];
                end
                
                % remove trials that dont match in block nr
                usetrls(data.blocknr(usetrls) ~= data.blocknr(trls)) = [];
                
                % fit logistic regression
                thisdat = data(usetrls, :);
                
                [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                    'binomial','link','logit');
                stimx = -4:0.01:4;
                yfit = glmval([b(1) - grandavg.overallLogistic(sj, 1); b(2)], stimx, 'logit');
                %  plot(x, yfit);
                
                % save betas
                grandavg.logistic(sj, r, u, whichLag, :) = b;
                grandavg.curve(sj, r, u, whichLag, :)    = yfit; % curve
                
            end % uncertainty bin
        end % resp
    end % sj
    
end

%% normalize by their overall bias
grandavg.logistic(:, :, :, :, whichLogisticParam) = bsxfun(@minus, grandavg.logistic(:, :, :, :, whichLogisticParam), grandavg.overallLogistic(:, whichLogisticParam));

% average over lags
lagGroups = {[1], [1:2], [1:3], [4:7]};

for l = 1:length(lagGroups),
    
    % take the mean over lags
    grandavg.groupedLogistic = squeeze(nanmean(grandavg.logistic(:, :, :, lagGroups{l}, :), 4)); % average over lags
    grandavg.groupedCurve = squeeze(nanmean(grandavg.curve(:, :, :, lagGroups{l}, :), 4)); % average over lags
    
    spcnt = spcnt + 1;
    
    %% plot 1. psychometric functions
    subplot(4,4,spcnt);     spcnt = spcnt + 1;
    
    hold on;
    
    % first plot the main psychometric function
    % plot(stimx, glmval(mean(grandavg.overallLogistic)', stimx, 'logit'), 'k');
    
    cnt = 1;
    colors = cbrewer('div', 'RdYlBu', 6);
    for ri = 1:2,
        for ui = 1,
            plot(stimx, squeeze(mean(grandavg.groupedCurve(:, ri, ui, :))), 'color', colors(cnt, :));
            cnt = cnt + 1;
        end
    end
    
    xlim([-3 3]); axis square;
    xlabel('Sensory evidence');
    ylabel('P_{resp}(A > B)');
    ylabel(sprintf('Lags %d, %d, %d, %d', lagGroups{l}));
    
    %% plot 2. plot interaction effect and 2way anova
    colors = linspecer(4);
    colors = colors([1 4], :);
    stimx2 = 1:u;
    
    subplot(4,4,spcnt); spcnt = spcnt + 1;
    hold on;
    plot([1 max(stimx2)], [0 0], '-', 'color', [0.4 0.4 0.4], 'LineWidth', 0.1);
    
    for r = [1 2],
        h = ploterr(stimx2, ...
            squeeze(nanmean(grandavg.groupedLogistic(:, r, :, whichLogisticParam))),  [], ...
            squeeze(nanstd(grandavg.groupedLogistic(:, r, :, whichLogisticParam))) ./ sqrt(length(subjects)), ...
            '-',  'hhxy', 0.001);
        set(h(1), 'color', colors(r, :), ...
            'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
        set(h(2), 'color', colors(r, :));
        stimx2 = stimx2 - 0.15;
    end
    
    switch whichLogisticParam
        case 1
            ylabel({'Response bias_{t0}'});
            ylim([-.2 .2]);
        case 2
            ylabel({'Logistic slope_{t0}'});
    end
    set(gca, 'ytick', [-.15 0  0.15]);
    whichModTitle = whichMod; whichModTitle(1) = upper(whichModTitle(1));
    axis square;
    xlim([0.5 3.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});
    % text(2, -.12, sprintf('Choice_{t-%d} B', whichLag), 'color', colors(1, :));
    % text(2, .12, sprintf('Choice_{t-%d} A', whichLag), 'color', colors(2, :));
    
    %% repeated measures anova on those bins
    cnt = 0; clear x s f
    for sj = 1:27,
        for r = 1:2,
            for u = 1:3,
                cnt = cnt + 1;
                x(cnt) = grandavg.groupedLogistic(sj, r, u, whichLogisticParam);
                s(cnt) = sj;
                f{1}(cnt) = r; % factor of previous response
                f{2}(cnt) = u; % factor of previous uncertainty
            end
        end
    end
    
    stats = rm_anova(x, s, f);
    
    %    >> stats = rm_anova(x, s, f);
    %           x = 1xN vector with data values
    %           s = 1xN vector with subject numbers, or the factor to repeat
    %           measures over
    %           f = factors, 1xNrFactors cell array. Each cell array should
    %           have an 1xN vector with the level within this factor.
    
    disp(stats.f1xf2)
    % title({sprintf('F_{(%d,%d)} = %.2f', stats.f1xf2.df(1), stats.f1xf2.df(2), stats.f1xf2.fstats); ...
    %    sprintf('p = %.3f', stats.f1xf2.pvalue)});
    
    %% plot 3. plot main effect on response repetition and one-way anova
    
    resp1 = -squeeze(grandavg.groupedLogistic(:, 1, :, whichLogisticParam));
    resp2 = squeeze(grandavg.groupedLogistic(:, 2, :, whichLogisticParam));
    grandavg.logisticRepGr = (resp1 + resp2) ./ 2;
    
    colors = linspecer(5);
    stimx2 = 1:u;
    
    subplot(4,4,spcnt); hold on; spcnt = spcnt + 1;
    h = ploterr(stimx2, ...
        squeeze(nanmean(grandavg.logisticRepGr)),  [], ...
        squeeze(nanstd(grandavg.logisticRepGr)) ./ sqrt(length(subjects)), ...
        '-',  'hhxy', 0.001);
    set(h(1), 'color', colors(5, :), ...
        'marker', '.', 'markerfacecolor', colors(5, :), 'markersize', 12);
    set(h(2), 'color', colors(5, :));
    plot([1 max(stimx2)], [0 0], '-', 'color', [0.4 0.4 0.4], 'LineWidth', 0.1);
    
    switch whichLogisticParam
        case 1
            ylabel({'Repetition bias_{t0}'});
            ylim([-.07 .15]);
        case 2
            ylabel({'Slope_{t0}'});
    end
    xlabel(sprintf('%s_{t-%d}',whichModTitle, whichLag));
    axis square;
    xlim([0.5 3.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});
    
    if 0,
        %% repeated measures anova on those bins
        cnt = 0; clear x s f
        for sj = 1:27,
            for u = 1:3,
                cnt = cnt + 1;
                x(cnt) = grandavg.logisticRepGr(sj, u);
                s(cnt) = sj;
                f{1}(cnt) = u; % factor of previous uncertainty
            end
        end
        clear stats
        stats = rm_anova(x, s, f);
        
        %    >> stats = rm_anova(x, s, f);
        %           x = 1xN vector with data values
        %           s = 1xN vector with subject numbers, or the factor to repeat
        %           measures over
        %           f = factors, 1xNrFactors cell array. Each cell array should
        %           have an 1xN vector with the level within this factor.
        
    end
    
    % disp(stats.f1)
    subplot(4,4,spcnt);
    text(0, 0.9, sprintf('Main F_{(%d,%d)} = %.2f', stats.f1.df(1), stats.f1.df(2), stats.f1.fstats));
    text(0, 0.75,  sprintf('p = %.3f', stats.f1.pvalue));
    
    text(0, 0.5, sprintf('Main modul F_{(%d,%d)} = %.2f', stats.f2.df(1), stats.f2.df(2), stats.f2.fstats));
    text(0, 0.35,  sprintf('p = %.3f', stats.f2.pvalue));
    
    text(0, 0.2, sprintf('Interaction F_{(%d,%d)} = %.2f', stats.f1xf2.df(1), stats.f1xf2.df(2), stats.f1xf2.fstats));
    text(0, 0.05,  sprintf('p = %.3f', stats.f1xf2.pvalue));
    axis off;
end

%% save fig
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/Fig4_valueShift_groupedLags_%s_%d_RTout%d.pdf', whichMod, whichLogisticParam, regressoutRT));
