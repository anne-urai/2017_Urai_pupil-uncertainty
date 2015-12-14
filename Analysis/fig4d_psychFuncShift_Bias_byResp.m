function fig4d_psychFuncShift_Bias_byResp(lagGroups, whichmodulator, grouping)

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end

% plot both the effect of pupil on overall repetition bias and show that
% this is symmetrical for both previous choices
subplot(441); hold on;

switch whichmodulator
    case 'pupil'
        whichMod = 'decision_pupil';
    case 'feedbackpupil'
        whichMod = 'feedback_pupil';
    case 'fb-decpupil'
        % see below for projecting out
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
        
        % in this case, take out decision effects
        switch whichmodulator
            case 'fb-decpupil'
                data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        end
        % outcome vector need to be 0 1 for logistic regression
        data.resp(data.resp == -1) = 0;
        
        % get an overall logistic fit for normalization
        grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, ...
            data.resp, ...
            'binomial','link','logit');
        
        % previous response
        resps = [0 1];
        for r = 1:2,
            
            % select for response and correctness
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
% plot for all subjects
% ========================================================= %
theseSj = 1:27;
colors = linspecer(4);
colors = colors([1 4], :);
stimx2 = 1:nbins;
plot([1 nbins], [0 0], 'k');

for r = [1 2],
    h = ploterr(stimx2, ...
        squeeze(nanmean(grandavg.logistic(theseSj, r, :, 1))),  [], ...
        squeeze(nanstd(grandavg.logistic(theseSj, r, :, 1))) ./ sqrt(length(theseSj)), ...
        '-',  'hhxy', 0.00001);
    set(h(1), 'color', colors(r, :), ...
        'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
    set(h(2), 'color', colors(r, :));
    stimx2 = stimx2 - 0.15;
    
    errorbar(nbins+1, nanmean(grandavg.respBiasNoPupilSplit(theseSj, r)), ...
        nanstd(grandavg.respBiasNoPupilSplit(theseSj, r)) ./ sqrt(length(theseSj)), ...
        '.', 'markersize', 12, 'color', colors(r, :));
end

switch whichMod,
    case 'baseline_pupil';
        xlabel(sprintf('Baseline pupil', 0));
    case 'decision_pupil'
        xlabel(sprintf('Pupil response'));
    case 'rt'
        xlabel('Reaction time');
    otherwise
        %  xlabel(sprintf('%s_{t-%d}', whichModTitle, whichLag));
end
axis tight;
xlim([0.5 nbins+1.5]); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high', 'all'});
ylabel('Next trial P(resp A)');

title('All subjects', 'color', colors(2, :));
% ========================================================= %
% plot for different subgroups
% ========================================================= %

load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'pupil'));

posRespSj = find(dat.response(:, 1) > 0);
negRespSj = find(dat.response(:, 1) < 0);

% ========================================================= %
% show that this effect is symmetrical
% ========================================================= %

for groups = 1:2,
    stimx2 = 1:nbins;
    switch groups
        case 1
            theseSj = posRespSj;
            subplot(442);
        case 2
            theseSj = negRespSj;
            subplot(443);
    end
    
    hold on;
    colors = linspecer(4);
    colors = colors([1 4], :);
    plot([1 nbins], [0 0], 'k');
    
    for r = [1 2],
        h = ploterr(stimx2, ...
            squeeze(nanmean(grandavg.logistic(theseSj, r, :, 1))),  [], ...
            squeeze(nanstd(grandavg.logistic(theseSj, r, :, 1))) ./ sqrt(length(theseSj)), ...
            '-',  'hhxy', 0.00001);
        set(h(1), 'color', colors(r, :), ...
            'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
        set(h(2), 'color', colors(r, :));
        stimx2 = stimx2 - 0.15;
        
        errorbar(nbins+1, nanmean(grandavg.respBiasNoPupilSplit(theseSj, r)), ...
            nanstd(grandavg.respBiasNoPupilSplit(theseSj, r)) ./ sqrt(length(theseSj)), ...
            '.', 'markersize', 12, 'color', colors(r, :));
    end
    
    switch whichMod,
        case 'baseline_pupil';
            xlabel(sprintf('Baseline pupil', 0));
        case 'decision_pupil'
            xlabel(sprintf('Pupil response'));
        case 'rt'
            xlabel(sprintf('Reaction time'));
        case 'feedback_pupil',
            xlabel('Feedback pupil');
            
            switch whichmodulator
                case 'fb-decpupil'
                    xlabel('Feedback-dec pupil');
            end
        otherwise
            %  xlabel(sprintf('%s_{t-%d}', whichModTitle, whichLag));
    end
    axis tight;
    xlim([0.5 nbins+1.5]); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high', 'all'});
    ylabel('Next trial P(resp A)');
    
    % repeated measures anova on those bins
    cnt = 0; clear x s f
    for sj = 1:length(theseSj),
        for u = 1:nbins,
            for r = 1:2,
                cnt = cnt + 1;
                x(cnt) = grandavg.logistic(theseSj(sj), r, u, 1, 1);
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
    % [~, pval(3)] = ttest(grandavg.respBiasNoPupilSplit(theseSj, 1), grandavg.respBiasNoPupilSplit(theseSj, 2));
    
    sigstar2({[1 nbins], [nbins+1 nbins+1]}, [stats.f1xf2.pvalue stats.f1.pvalue]);
    stats.f1xf2
    
    colors = linspecer(9);
    
    switch groups
        case 1
            title('Repeaters', 'color', colors(8,:));
        case 2
            title('Switchers', 'color', colors(9,:));
            
            % add text to indicate whats going on
            colors = linspecer(4);
            text(nbins+1.5, 0.03, 'Resp_{-1} A', 'color', colors(4, :));
            text(nbins+1.5, -0.03, 'Resp_{-1} B', 'color', colors(1,:));
    end
end

% suplabel(sprintf('Lags %d %d %d', lagGroups), 'x');
% print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4c_biasShift_%s_correct%d_lags%d%d%d.pdf', whichMod, correctness, lagGroups));

end
