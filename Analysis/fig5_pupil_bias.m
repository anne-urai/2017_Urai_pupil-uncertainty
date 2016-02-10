function [grandavg] = fig5_pupil_bias(nbins)
% plots uncertainty by accuracy both for the modelfits and the pupil

global mypath;
nbins = 3;
cnt = 1;
%% fit logistic slope on high vs low pupil bins

pupfields = {'earlydecision_pupil', 'decision_pupil'};
for p = 1:length(pupfields),
    
    subjects = 1:27;
    grandavg.betas = nan(length(subjects), nbins, 2);
    
    for sj = unique(subjects),
        % get all the data
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        
        pupQs = quantile(data.(pupfields{p}), nbins - 1);
        % divide into low and high pupil bins
        % puptrls{1} = find(data.(pupfields{p}) < quantile(data.(pupfields{p}), 0.5));
        % puptrls{2} = find(data.(pupfields{p}) > quantile(data.(pupfields{p}), 0.5));
        
        for b = 1:nbins,
            
            switch b
                case 1
                    trls = find(data.(pupfields{p}) < pupQs(b));
                case nbins
                    trls = find(data.(pupfields{p}) > pupQs(b-1));
                otherwise
                    trls = find(data.(pupfields{p}) > pupQs(b-1) & data.(pupfields{p}) < pupQs(b));
            end
            
            laggedtrls = trls + 1;
            
            % exclude trials at the end of the block
            if any(laggedtrls > size(data, 1)),
                trls(laggedtrls > size(data, 1)) = [];
                laggedtrls(laggedtrls > size(data, 1)) = [];
            end
            
            % remove trials that dont match in block nr
            removeTrls = data.blocknr(laggedtrls) ~= data.blocknr(trls);
            laggedtrls(removeTrls) = [];
            trls(removeTrls) = [];
            
            
            thisdat = data(laggedtrls, :);
            
            % fit on all the trials, dont average beforehand!
            x = thisdat.motionstrength;
            y = thisdat.resp; y(y==-1) = 0;
            
            % sort both
            [x, idx] = sort(x);
            y = y(idx);
            
            [beta,dev,stats] = glmfit(x, [y ones(size(y))], ...
                'binomial','link','logit');
            % newx = -3:0.01:3;
            % yfit = glmval(beta, newx,'logit');
            % plot(newx, yfit,'-', 'color', cols(b, :));
            
            % also save the curve itself
            grandavg.betas(sj, b, :) = beta;
            % grandavg.yfit(sj, b, :) = yfit;
            
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
    subplot(4,4,cnt); cnt = cnt + 1;
    hold on;
    bar(1, mean(grandavg.betas(:, 1, 2)), 'facecolor', [0.6 0.6 0.6], 'edgecolor', 'none', 'barwidth', 0.4);
    bar(2, mean(grandavg.betas(:, 2, 2)), 'facecolor', [0.4 0.4 0.4], 'edgecolor', 'none', 'barwidth', 0.4);
    bar(3, mean(grandavg.betas(:, 3, 2)), 'facecolor', [0.2 0.2 0.3], 'edgecolor', 'none', 'barwidth', 0.4);
    
    errorbar(1:nbins, squeeze(nanmean(grandavg.betas(:, :, 2))),  ...
        squeeze(std(grandavg.betas(:, :, 2))) / sqrt(length(subjects)), '.k');
    
    % do ttest on regression coefficients
    [~, pval(3), ~, stat] = ttest(grandavg.betas(:, 1, 2), grandavg.betas(:, end, 2));
    xlim([0.5 nbins+0.5]); box off;
    mysigstar([1 nbins], [1 1], pval(3));
    ylim([0.4 1.1]);
    set(gca, 'ytick', [0.5 1]);
    xlabel('Pupil response'); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high'});
    ylabel('Next trial slope');
    title(pupfields{p}, 'interpreter', 'none');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    
    %% absolute bias as a function of pupil!
    subplot(4,4,cnt); cnt = cnt + 1;
    hold on;
    absb = abs(squeeze(grandavg.betas(:, :, 1)));
    bar(1, mean(absb(:, 1)), 'facecolor', [0.6 0.6 0.6], 'edgecolor', 'none', 'barwidth', 0.4);
    bar(2, mean(absb(:, 2)), 'facecolor', [0.4 0.4 0.4], 'edgecolor', 'none', 'barwidth', 0.4);
    bar(3, mean(absb(:, 3)), 'facecolor', [0.2 0.2 0.2], 'edgecolor', 'none', 'barwidth', 0.4);
    
    errorbar(1:3, squeeze(nanmean(absb)),  ...
        squeeze(std(absb)) / sqrt(length(subjects)), '.k');
    
    % do ttest on regression coefficients
    [~, pval(3), ~, stat] = ttest(absb(:, 1), absb(:, end));
    xlim([0.5 nbins+0.5]); box off;
    mysigstar([1 nbins], [0.45 0.45], pval(3));
    ylim([0.1 0.5]);
    xlabel('Pupil response'); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high'});
    ylabel('Next bias size');
    title(pupfields{p}, 'interpreter', 'none');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    
end

end
