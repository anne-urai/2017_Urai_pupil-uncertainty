function fig4d_psychFuncShift_Bias(lagGroups)

if ~exist('lagGroups', 'var'), lagGroups = 1; end

% plot both the effect of pupil on overall repetition bias and show that
% this is symmetrical for both previous choices
subplot(445);

whichmodulator = 'pupil';

switch whichmodulator
    case 'pupil'
        whichMod = 'feedback_pupil';
    case 'rt'
        whichMod = 'rt';
end

nbins = 3;
subjects = 1:27;
clear grandavg;
whichLags = 1:3; % lag 1

for lag = whichLags,
    for sj = unique(subjects),
        data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
        switch whichMod
            case 'rt'
                data.rt = log(data.rt);
            case 'decision_pupil';
                %   data.decision_pupil = projectout(data.decision_pupil, data.rt);
        end
        
        % outcome vector need to be 0 1 for logistic regression
        data.resp(data.resp == -1) = 0;
        
        % get an overall logistic fit
        grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, data.resp, ...
            'binomial','link','logit');
        
        % previous response
        resps = [0 1];
        for r = 1:2,
            
            trls = find(data.resp == resps(r));
            laggedtrls = trls+lag;
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
            
            % get an overall logistic fit
            grandavg.respBiasNoPupilSplit(sj, r, lag, :) = glmfit(thisdat.motionstrength, thisdat.resp, ...
                'binomial','link','logit');
            
            %% split into pupil bins
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
                laggedtrls = trls+lag;
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
                grandavg.logistic(sj, r, u, lag, :) = b;
                
            end % uncertainty bin
        end % resp
    end % sj
end

% ========================================================= %
% normalize by their overall bias
% ========================================================= %

grandavg.logistic = bsxfun(@minus, grandavg.logistic, grandavg.overallLogistic(:, 1));
grandavg.respBiasNoPupilSplit = bsxfun(@minus, grandavg.respBiasNoPupilSplit, grandavg.overallLogistic(:, 1));

% ========================================================= %
% make one plot with repetition (rather than response) bias
% ========================================================= %

resp1 = squeeze(grandavg.logistic(:, 1, :, :, 1));
resp2 = squeeze(grandavg.logistic(:, 2, :, :, 1));

grandavg.logisticRep = (-resp1 + resp2) ./ 2;

resp1 = squeeze(grandavg.respBiasNoPupilSplit(:, 1, :, 1));
resp2 = squeeze(grandavg.respBiasNoPupilSplit(:, 2, :, 1));

grandavg.repetitionBias = (-resp1 + resp2) ./ 2;

% ========================================================= %
% combine the lag groups
% ========================================================= %

grandavg.logisticRep = squeeze(nanmean(grandavg.logisticRep(:, :, lagGroups), 3));
grandavg.logistic = squeeze(nanmean(grandavg.logistic(:, :, :, lagGroups, :), 4));
grandavg.respBiasNoPupilSplit = squeeze(nanmean(grandavg.respBiasNoPupilSplit(:, :, lagGroups, 1), 3));
grandavg.repetitionBias  = squeeze(nanmean(grandavg.repetitionBias(:, lagGroups), 2));

% ========================================================= %
% plot
% ========================================================= %

stimx2 = 1:u;
hold on;
h = ploterr(stimx2, ...
    nanmean(grandavg.logisticRep),  [], ...
    nanstd(grandavg.logisticRep) ./ sqrt(length(subjects)), ...
    '-k',  'hhxy', 0.001);
set(h(1), 'color', 'k', ...
    'marker', '.', 'markerfacecolor', 'k', 'markersize', 12);
set(h(2), 'color', 'k');

% next to it, plot the overall bias
errorbar(4, nanmean(grandavg.repetitionBias), ...
    nanstd(grandavg.repetitionBias) ./ sqrt(length(subjects)), ...
    '.k', 'markersize', 12);
axis square; axis tight;
xlim([0.5 nbins+1.5]); set(gca, 'xtick', 1:4, ...
    'xticklabel', {'low', 'med', 'high', 'all'}, 'xticklabelrotation', 0);
box off;

switch whichMod,
    case 'baseline_pupil';
        xlabel(sprintf('Baseline pupil_{t%d}', 0));
    case 'decision_pupil'
        xlabel(sprintf('Decision pupil_{-1}'));
    otherwise
      %  xlabel(sprintf('%s_{t-%d}', whichModTitle, whichLag));
end
ylabel('Repetition bias_{0}');


% repeated measures anova on those bins
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
[~, pval(1)] = ttest(grandavg.logisticRep(:, 1), grandavg.logisticRep(:, 2));
[~, pval(2)] = ttest(grandavg.logisticRep(:, 2), grandavg.logisticRep(:, 3));

[pval(3)] = randtest1d(grandavg.repetitionBias);
sigstar({[1 3], [4 4]}, [stats.f1.pvalue pval(3)]);

% ========================================================= %
% show that this effect is symmetrical
% ========================================================= %

subplot(446); hold on;
colors = linspecer(4);
colors = colors([1 4], :);
plot([1 3], [0 0], 'k');
for r = [1 2],
    h = ploterr(stimx2, ...
        squeeze(nanmean(grandavg.logistic(:, r, :, 1))),  [], ...
        squeeze(nanstd(grandavg.logistic(:, r, :, 1))) ./ sqrt(length(subjects)), ...
        '-',  'hhxy', 0.00001);
    set(h(1), 'color', colors(r, :), ...
        'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
    set(h(2), 'color', colors(r, :));
    stimx2 = stimx2 - 0.15;
    
    errorbar(4, nanmean(grandavg.respBiasNoPupilSplit(:, r)), ...
        nanstd(grandavg.respBiasNoPupilSplit(:, r)) ./ sqrt(length(subjects)), ...
        '.', 'markersize', 12, 'color', colors(r, :));
end

switch whichMod,
    case 'baseline_pupil';
        xlabel(sprintf('Baseline pupil_{t%d}', 0));
    case 'decision_pupil'
        xlabel(sprintf('Decision pupil_{-1}'));
    otherwise
      %  xlabel(sprintf('%s_{t-%d}', whichModTitle, whichLag));
end
axis tight;
axis square;
xlim([0.5 nbins+1.5]); set(gca, 'xtick', 1:4, 'xticklabel', {'low', 'med', 'high', 'all'});
ylabel('Response bias_{0}');

text(2, 0.1, 'Resp_{-1} A', 'color', colors(2, :));
text(2, -0.1, 'Resp_{-1} B', 'color', colors(1,:));

% repeated measures anova on those bins
cnt = 0; clear x s f
for sj = 1:27,
    for u = 1:3,
        for r = 1:2,
            cnt = cnt + 1;
            x(cnt) = grandavg.logistic(sj, r, u, 1, 1);
            s(cnt) = sj;
            f{1}(cnt) = r; % factor of previous response
            f{2}(cnt) = u; % factor of previous uncertainty
        end
    end
end
clear stats
stats = rm_anova(x, s, f);

% also test for the 2 pairwise comparisons
% [h, pval(1)] = ttest(grandavg.logistic(:, 1, whichParam), grandavg.logistic(:, 2, whichParam));
% [h, pval(2)] = ttest(grandavg.logistic(:, 2, whichParam), grandavg.logistic(:, 3, whichParam));

[pval(3)] = randtest1d(grandavg.respBiasNoPupilSplit(:, 1), grandavg.respBiasNoPupilSplit(:, 2));
%sigstar({[1 3], [4 4]}, [stats.f1xf2.pvalue pval(3)]);

suplabel(sprintf('Lags %d %d %d', lagGroups), 'x');
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4c_psychFuncShift_bias_lags%d%d%d.pdf', lagGroups));

end
