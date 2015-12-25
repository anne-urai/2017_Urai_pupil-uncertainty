function fig4d_psychFuncShift_Slope
% plot both the effect of pupil on overall repetition bias and show that
% this is symmetrical for both previous choices

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
whichLag = 1; % lag 1

for sj = unique(subjects),
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    switch whichMod
        case 'rt'
            data.rt = log(data.rt);
        case 'decision_pupil';
            data.decision_pupil = projectout(data.decision_pupil, data.rt);
            
    end
    
    % outcome vector need to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    
    % get an overall logistic fit
    grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, data.resp, ...
        'binomial','link','logit');
    
    % previous response
    if nbins > 2,
        uncQs = quantile(data.(whichMod), nbins - 1);
    elseif nbins == 2,
        uncQs = median(data.(whichMod));
    end
    
    % uncertainty bins
    for u = 1:nbins,
        
        switch u
            case 1
                trls = find(data.(whichMod) < uncQs(u));
            case nbins
                trls = find(data.(whichMod) > uncQs(u-1));
            otherwise
                trls = find(data.(whichMod) > uncQs(u-1) & data.(whichMod) < uncQs(u));
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
        
        [b, d] = glmfit(thisdat.motionstrength, thisdat.resp, ...
            'binomial','link','logit');
        
        % save betas
        grandavg.logistic(sj, u, whichLag, :) = b;
        
    end % uncertainty bin
end % sj

% ========================================================= %
% normalize by their overall slope
% ========================================================= %

grandavg.logistic = bsxfun(@minus, grandavg.logistic, grandavg.overallLogistic(:, 2));

% ========================================================= %
% one plot with slope
% ========================================================= %

grandavg.logisticSlope = squeeze(grandavg.logistic(:, :, 1, 2));

stimx2 = 1:u;
subplot(4,4,1);
h = ploterr(stimx2, ...
    nanmean(grandavg.logisticSlope),  [], ...
    nanstd(grandavg.logisticSlope) ./ sqrt(length(subjects)), ...
    '-k',  'hhxy', 0.001);
set(h(1), 'color', 'k', ...
    'marker', '.', 'markerfacecolor', 'k', 'markersize', 12);
set(h(2), 'color', 'k');

switch whichMod,
    case 'baseline_pupil';
        xlabel(sprintf('Baseline pupil_{t%d}', 0));
    case 'decision_pupil'
        xlabel(sprintf('Decision pupil_{t-1}'));
    otherwise
        xlabel(sprintf('%s_{t-%d}', whichModTitle, whichLag));
end
ylabel({'Slope_{t0}'; '(change from avg'});

axis square;
xlim([0.5 3.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});
ylim([-0.02 0.04]); box off;

% repeated measures anova on those bins
cnt = 0; clear x s f
for sj = 1:27,
    for u = 1:3,
        cnt = cnt + 1;
        x(cnt) = grandavg.logisticSlope(sj, u);
        s(cnt) = sj;
        f{1}(cnt) = u; % factor of previous uncertainty
    end
end
clear stats
stats = rm_anova(x, s, f);

% also test for the 2 pairwise comparisons
[~, pval(1)] = ttest(grandavg.logisticSlope(:, 1), grandavg.logisticSlope(:, 2));
[~, pval(2)] = ttest(grandavg.logisticSlope(:, 2), grandavg.logisticSlope(:, 3));

sigstar({[1 3]}, [stats.f1.pvalue]);

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4d_psychFuncShift_slope.pdf'));

end
