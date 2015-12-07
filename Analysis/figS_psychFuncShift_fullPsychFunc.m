function [] = fig4c_psychFuncShift_interaction()
% ========================================================= %
% panel C: shift of the mean of the psychometric function on lag 1-3 with complete stats anova
% ========================================================= %

clf; clear; clc;
whichmodulator = 'pupil';

switch whichmodulator
    case 'pupil'
        whichMod = 'decision_pupil';
    case 'rt'
        whichMod = 'rt';
end

nbins = 3;
subjects = 1:27;
clear grandavg;
lagGroups = 1; % lag 1
spcnt = 1;

for whichLag = 1,
    for sj = unique(subjects),
        data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
        switch whichMod
            case 'rt'
                data.rt = log(data.rt);
        end
        
        % outcome vector need to be 0 1 for logistic regression
        data.resp(data.resp == -1) = 0;
        % remove the trials at the end of each block
        endOfBlockTrls = find(data.trialnr == 50);
        
        % get an overall logistic fit
        b = glmfit(data.motionstrength, data.resp, ...
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
                laggedtrls = trls+whichLag;
                % exclude trials at the end of the block
                if any(laggedtrls > size(data, 1)),
                    trls(laggedtrls > size(data, 1)) = [];
                    laggedtrls(laggedtrls > size(data, 1)) = [];
                end
                
                % remove trials that dont match in block nr
                removeTrls = data.blocknr(laggedtrls) ~= data.blocknr(trls);
                laggedtrls(removeTrls) = [];
                trls(removeTrls) = [];
                
                % fit logistic regression
                thisdat = data(laggedtrls, :);
                
                [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                    'binomial','link','logit');
                
                % save betas
                grandavg.logistic(sj, r, u, whichLag, :) = b;
                
                % also fit a function with lapse rate included
                logisticFit = Psychometric_Logistic(thisdat, 0,0);
                grandavg.logisticLapse(sj, r, u, whichLag, :) = logisticFit.mean.beta;
                
            end % uncertainty bin
        end % resp
    end % sj
end

for whichParam = 1:4,
    
    % get the parameter we want
    grandavg.logistic = grandavg.logisticLapse(:, :, :, :, whichParam);
    
    % take the mean over lags
    grandavg.groupedLogistic = squeeze(nanmean(grandavg.logistic(:, :, :, lagGroups, :), 4)); % average over lags
    
    % compute main effect
    resp1 = -squeeze(grandavg.logistic(:, 1, :, lagGroups));
    resp2 = squeeze(grandavg.logistic(:, 2, :, lagGroups));
    
    grandavg.logisticRep = (resp1 + resp2) ./ 2;
    stimx2 = 1:u;
    
    subplot(4,4,spcnt); spcnt = spcnt + 1;
    hold on;
    h = ploterr(stimx2, ...
        squeeze(nanmean(grandavg.logisticRep)),  [], ...
        squeeze(nanstd(grandavg.logisticRep)) ./ sqrt(length(subjects)), ...
        '-',  'hhxy', 0.001);
    set(h(1), 'color', 'k', ...
        'marker', '.', 'markerfacecolor', 'k', 'markersize', 12);
    set(h(2), 'color', 'k');
   % plot([1 max(stimx2)], [0 0], '-', 'color', [0.4 0.4 0.4], 'LineWidth', 0.1);
    
    switch whichParam
        case 1
            ylabel({'Repetition bias_{t0}'});
        case 2
            ylabel({'Slope_{t0}'});
        case 3
            ylabel({'GuessRate_{t0}'});
        case 4
            ylabel({'LapseRate_{t0}'});
    end
    
   % ylim([-.07 .15]);
    switch whichMod,
        case 'baseline_pupil';
            xlabel(sprintf('Baseline pupil_{t%d}', 0));
        case 'decision_pupil'
            xlabel(sprintf('Decision pupil_{t-1}'));
        otherwise
            xlabel(sprintf('%s_{t-%d}', whichModTitle, whichLag));
    end
    axis square;
    xlim([0.5 3.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});
    
    %% repeated measures anova on those bins
    cnt = 0; clear x s f
    for sj = 1:27,
        for u = 1:3,
            cnt = cnt + 1;
            x(cnt) = grandavg.logisticRep(sj, u);
            s(cnt) = sj;
            f{1}(cnt) = u; % factor of previous uncertainty
        end
    end
    clear stats
    stats = rm_anova(x, s, f);
    
    % also test for the 2 pairwise comparisons
    [h, pval(1)] = ttest(grandavg.logisticRep(:, 1), grandavg.logisticRep(:, 2));
    [h, pval(2)] = ttest(grandavg.logisticRep(:, 2), grandavg.logisticRep(:, 3));
    
    sigstar({[1 2], [2 3], [1 3]}, [pval(1) pval(2) stats.f1.pvalue]);
    
    % ========================================================= %
    %% now show the full interaction
    % ========================================================= %
    colors = linspecer(9, 'qualitative');
    colors = colors([5 9], :); % grey and purple
    stimx2 = 1:u;
    
    subplot(4,4,spcnt); spcnt = spcnt +3;
    hold on;
   % plot([1 max(stimx2)], [0 0], '-', 'color', 'k', 'LineWidth', 0.1);
    
    for r = [1 2],
        h = ploterr(stimx2, ...
            squeeze(nanmean(grandavg.groupedLogistic(:, r, :, 1))),  [], ...
            squeeze(nanstd(grandavg.groupedLogistic(:, r, :, 1))) ./ sqrt(length(subjects)), ...
            '-',  'hhxy', 0.00001);
        set(h(1), 'color', colors(r, :), ...
            'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
        set(h(2), 'color', colors(r, :));
        stimx2 = stimx2 - 0.15;
    end
    
    switch whichParam
        case 1
            ylabel({'Bias_{t0}'});
        case 2
            ylabel({'Slope_{t0}'});
        case 3
            ylabel({'GuessRate_{t0}'});
        case 4
            ylabel({'LapseRate_{t0}'});
    end
    
    set(gca, 'ytick', [-.15 0  0.15]);
    xlim([0.5 3.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});
    % ylim([-.15 .15]);
    switch whichmodulator
        case 'rt'
            xlabel('Reaction time_{t-1}');
        otherwise
            xlabel('Decision pupil_{t-1}');
    end
    axis square;
    
    % repeated measures anova on those bins
    cnt = 0; clear x s f
    for sj = 1:27,
        for r = 1:2,
            for u = 1:3,
                cnt = cnt + 1;
                x(cnt) = grandavg.groupedLogistic(sj, r, u, 1);
                s(cnt) = sj;
                f{1}(cnt) = r; % factor of previous response
                f{2}(cnt) = u; % factor of previous uncertainty
            end
        end
    end
    
    stats = rm_anova(x, s, f);
    sigstar({[1 3]}, stats.f1xf2.pvalue);
    
    % also test for the 2 pairwise comparisons
    % [h, pval(1)] = ttest(grandavg.logistic(:, 1, 1), grandavg.logisticRep(:, 2, 1));
    % [h, pval(2)] = ttest(grandavg.logisticRep(:, 2), grandavg.logisticRep(:, 3));
    
    % sigstar({[1 2], [2 3], [1 3]}, [pval(1) pval(2) stats.f1xf2.pvalue]);
end


print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4c_psychFuncShift_FullPsychFunc.pdf'));



