function [grandavg] = fig5_pupil_bias()
% plots uncertainty by accuracy both for the modelfits and the pupil

global mypath;
nbins = 3;
cnt = 1;
%% fit logistic slope on high vs low pupil bins

pupfields = {'rt'};
for p = 1:length(pupfields),
    
    subjects = 1:27;
    grandavg.betas = nan(length(subjects), nbins, 2);
    
    for sj = unique(subjects),
        % get all the data
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        
        % find their individual 70% correct threshold
        xrange          = abs(data.motionstrength);
        yrange          = data.correct;
        [newx, newy]    = divideintobins(xrange, yrange, 10);
        % cumulative normal at 0.5
        ft              = fittype( @(b, s, x) 0.5 * (1 + erf( (abs(x-b)/(s * sqrt(2))))));
        fitobj          = fit(xrange, yrange, ft);
        % fit this to these data
        sortedx         = sort(xrange);
        yval            = feval(fitobj, sortedx);
        threshold       = sortedx(dsearchn(yval, 0.6));
        
        pupQs = quantile(data.(pupfields{p}), nbins - 1);
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
            grandavg.betas(sj, b, :) = beta;
            
            % save the accuracy
            grandavg.accuracy(sj, b) = nanmean(thisdat.correct(thisdat.motionstrength < threshold));
        end
    end
    
    % make sure I actually plot the bar graph here to reproduce figure 3!
    subplot(4,4,cnt); cnt = cnt + 1;
    errorbar(1:nbins, squeeze(nanmean(grandavg.accuracy)),  ...
        squeeze(std(grandavg.accuracy)) / sqrt(length(subjects)), '.k-', 'markersize', 12);
    
    % do ttest on regression coefficients
    [~, pval(3), ~, stat] = ttest(grandavg.accuracy(:, 1), grandavg.accuracy(:, end));
    xlim([0.5 nbins+0.5]); box off;
    mysigstar([1 nbins], [0.74 0.74], pval(3));
    set(gca, 'ytick', [0.68 0.7 0.72 0.74]);
    ylim([0.68 0.74]);
    xlabel('Pupil response'); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high'});
    ylabel('Next trial P(correct)');
    
    %%  bias as a function of pupil!
    subplot(4,4,cnt); cnt = cnt + 1;
    absb = squeeze(grandavg.betas(:, :, 1));
    errorbar(1:nbins, squeeze(nanmean(absb)),  ...
        squeeze(std(absb)) / sqrt(length(subjects)), '.k-', 'markersize', 12);
    % do ttest on regression coefficients
    [~, pval(3), ~, stat] = ttest(absb(:, 1), absb(:, end));
    xlim([0.5 nbins+0.5]); box off;
    bf10 = t1smpbf(stat.tstat,27)
    disp(1/bf10)
    mysigstar([1 nbins], [0.1 0.1], pval(3));
    ylim([-0.12 0.12]); set(gca, 'ytick', [-0.1 0 0.1]);
    xlabel('Pupil response'); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high'});
    ylabel('Next trial bias');
    
    %% absolute bias as a function of pupil!
    subplot(4,4,cnt); cnt = cnt + 1;
    absb = abs(squeeze(grandavg.betas(:, :, 1)));
    errorbar(1:nbins, squeeze(nanmean(absb)),  ...
        squeeze(std(absb)) / sqrt(length(subjects)), '.k-', 'markersize', 12);
    % do ttest on regression coefficients
    [~, pval(3), ~, stat] = ttest(absb(:, 1), absb(:, end));
    bf10 = t1smpbf(stat.tstat,27)
    disp(1/bf10)
    xlim([0.5 nbins+0.5]); box off;
    mysigstar([1 nbins], [0.4 0.4], pval(3));
    ylim([0.2 0.41]); set(gca, 'ytick', [0.2 0.3 0.4]);
    xlabel('Pupil response'); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high'});
    ylabel('Next trial |bias|');
    
end

end
