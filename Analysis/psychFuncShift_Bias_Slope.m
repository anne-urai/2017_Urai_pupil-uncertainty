function psychFuncShift_Bias_Slope(whichmodulator, nbins, correctness)
global mypath;

if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('correctness', 'var'); correctness = []; end
if ~exist('nbins', 'var'); nbins = 3; end
lag = 1; % look at 1 trial in the past

warning('error', 'stats:glmfit:PerfectSeparation');

switch whichmodulator
    case 'pupil'
        whichMod = 'decision_pupil';
    case 'fbpupil'
        whichMod = 'feedback_pupil';
    case 'fb+decpupil'
        % see below for projecting out
        whichMod = 'feedback_pupil';
    case 'dec+fbpupil'
        % see below for projecting out
        whichMod = 'decision_pupil';
    otherwise
        whichMod = 'uncertainty';
end

% =========================================== %
% BIAS - ITERATE OVER TWO RESPONSES
% =========================================== %

subjects = 1:27;
clear grandavg;

for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data((data.sessionnr > 1), :); % get rid of residual learning effects
    
    switch whichmodulator
        case 'fb+decpupil'
            % in this case, take out decision effects
            data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        case 'dec+fbpupil'
            data.decision_pupil = projectout(data.decision_pupil, data.feedback_pupil);
            
        case 'pupil',
            data.decision_pupil = projectout(data.decision_pupil, zscore(log(data.rt + 0.1)));
        case 'rt'
            data.rt = projectout(zscore(log(data.rt+0.1)), data.decision_pupil);
        case 'evidence'
            % single-trial evidence strength is absolute motionenergy
            data.evidence = abs(data.motionstrength);
            % make sure to use the same range of evidence for each session
            for session = unique(data.sessionnr)',
                data.evidence(data.sessionnr == session) = ...
                    zscore(data.evidence(data.sessionnr == session));
            end
    end
    
    % outcome vector need to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    
    % get an overall logistic fit for normalization
    grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, ...
        data.resp, 'binomial','link','logit');
    
    % previous response
    resps = [0 1];
    for r = 1:2,
        
        % split into quantiles for each response
        if isempty(correctness),
            uncQs = quantile(data.(whichMod)(data.resp == resps(r)), nbins - 1);
            % uncQs = median(data.(whichMod)(data.resp == resps(r)));
        else
            uncQs = quantile(data.(whichMod)(data.resp == resps(r) & data.correct == correctness), nbins - 1);
            % uncQs = median(data.(whichMod)(data.resp == resps(r) & data.correct == correctness));
        end
        
        % uncertainty bins
        for u = 1:nbins,
            
            % find the trials in this bin
            if isempty(correctness), % take all trials
                switch u
                    case 1
                        trls = find(data.resp == resps(r) & data.(whichMod) <= uncQs(u));
                    case nbins % last one
                        trls = find(data.resp == resps(r) & data.(whichMod) > uncQs(u-1));
                    otherwise
                        trls = find(data.resp == resps(r) & ...
                            data.(whichMod) > uncQs(u-1) & data.(whichMod) <= uncQs(u));
                end
            else % either correct or error trials
                switch u
                    case 1
                        trls = find(data.resp == resps(r) & data.correct == correctness & data.(whichMod) <= uncQs(u));
                    case nbins % last one
                        trls = find(data.resp == resps(r) & data.correct == correctness & data.(whichMod) > uncQs(u-1));
                    otherwise
                        trls = find(data.resp == resps(r) & data.correct == correctness & ...
                            data.(whichMod) > uncQs(u-1) & data.(whichMod) <= uncQs(u));
                end
            end
            
            % with this selection, take the trials after that
            laggedtrls = trls+lag;
            
            % exclude trial at the end of the dataset
            if any(laggedtrls > size(data, 1)),
                trls(laggedtrls > size(data, 1)) = [];
                laggedtrls(laggedtrls > size(data, 1)) = [];
            end
            
            % remove trials that dont match in block nr
            removeTrls = find(data.blocknr(laggedtrls) ~= data.blocknr(trls));
            laggedtrls(removeTrls) = [];
            
            % fit logistic regression
            thisdat = data(laggedtrls, :);
            
            try
                b = glmfit(thisdat.motionstrength, thisdat.resp, ...
                    'binomial','link','logit');
            catch
                warning('putting nan in betas!');
                b = nan(size(b));
            end
            
            % save betas
            grandavg.logistic(sj, r, u, :) = b;
            
        end % uncertainty bin
    end % resp
end % sj

% ========================================================= %
% normalize by their overall bias for A or B
% keep only the bias, discard slope info
% ========================================================= %

% grandavg.logistic = bsxfun(@minus, grandavg.logistic(:, :, :, 1), grandavg.overallLogistic(:, 1));

grandavg.logistic = grandavg.logistic(:, :, :, 1);

% ========================================================= %
% compute repetition (rather than response) bias
% ========================================================= %

resp1 = -1 * (squeeze(grandavg.logistic(:, 1, :)) - 1) + 1;
resp1 = -1 * squeeze(grandavg.logistic(:, 1, :));

resp2 = squeeze(grandavg.logistic(:, 2, :));

% since this is centred at 0.5, treat it that way
grandavg.logisticRep = (resp1 + resp2) ./ 2;

% ========================================================= %
% transform from log(odds) to probability
% ========================================================= %

grandavg.logisticRep = exp(grandavg.logisticRep)./(1+exp(grandavg.logisticRep));
grandavgBias = grandavg;

% =========================================== %
% SLOPE - DO NOT ITERATE OVER TWO RESPONSES   %
% =========================================== %

clearvars -except grandavgBias subjects mypath whichmodulator whichMod nbins correctness lag

for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data(find(data.sessionnr > 1), :);
    
    % in this case, take out decision effects
    
    switch whichmodulator
        case 'fb+decpupil'
            data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        case 'dec+fbpupil'
            data.decision_pupil = projectout(data.decision_pupil, data.feedback_pupil);
            
        case 'pupil',
            data.decision_pupil = projectout(data.decision_pupil, zscore(log(data.rt + 0.1)));
        case 'rt'
            data.rt = projectout(zscore(log(data.rt+0.1)), data.decision_pupil);
        case 'evidence'
            data.evidence = abs(data.motionstrength);
        case 'uncertainty';
            % fit a probit slope so we can get the sigma
            b       = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');
            sigma   = 1/b(2);  % standard deviation at these values is the inverse!
            bound   = -b(1); % the bound is negative if people say
            
            % for each trial, compute the average level of uncertainty
            data.uncertainty = arrayfun(@simulateUncertainty, abs(data.motionstrength), ...
                data.correct, sigma*ones(length(data.correct), 1), bound*ones(length(data.correct), 1));
    end
    
    % outcome vector needs to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    
    % previous response
    if ~isempty(correctness),
        uncQs = quantile(data.(whichMod)(data.correct == correctness), nbins - 1);
    else
        uncQs = quantile(data.(whichMod), nbins - 1);
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
        
        if ~isempty(correctness),
            % continue with a subset of trials
            trls = intersect(trls, find(data.correct == correctness));
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
        
        % fit logistic regression
        thisdat = data(laggedtrls, :);
        assert(~isempty(thisdat));
        grandavg.logistic(sj, u, :) = glmfit(thisdat.motionstrength, thisdat.resp, ...
            'binomial','link','logit');
        
    end % uncertainty bin
end % sj

% replace
grandavg.logisticSlope = squeeze(grandavg.logistic(:, :, 2));
grandavgSlope = grandavg;

% ========================================================= %
% plot
% ========================================================= %
%
x   = 1:nbins;
y1  = grandavgBias.logisticRep;
y2  = grandavgSlope.logisticSlope;

% Store the axes handles produced by PLOTYY
[ax, p1, p2] = plotyy(x, nanmean(y1), x, nanmean(y2));
p1.Color     = 'k';
p2.Color     = [0.4 0.4 0.4];

% Use HOLD and ERRORBAR, passing axes handles to the functions.
hold(ax(1), 'on');
errorbar(ax(1), x, nanmean(y1), nanstd(y1) ./sqrt(27), ...
    '.-', 'color', 'k', 'markerfacecolor', 'k', 'markersize', 12);

set(ax(1), 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins, ...
    'xticklabel', {'low', 'med', 'high'}, 'ylim', [0.48 0.56], 'ytick', [0.49 0.51 0.53], ...
    'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off');
% low vs high
[~, pval, ci, stats] = ttest(y1(:, 1), y1(:, end));
bf10 = t1smpbf(stats.tstat, 27);
fprintf('bias, bf10 = %f', bf10);
mysigstar(gca, [1 nbins], [0.53 0.53], pval, 'k', 'down');
axis(ax(1), 'square');

% second errorbar, grey
hold(ax(2), 'on');
errorbar(ax(2), x, nanmean(y2), nanstd(y2) ./sqrt(27), ...
    '.-', 'color', [0.4 0.4 0.4], 'markerfacecolor', [0.4 0.4 0.4], 'markersize', 12);

set(ax(2), 'xlim', [0.5 nbins+0.5], 'ylim', [0.40 0.9], 'ytick', [0.8 0.9], ...
    'xcolor', 'k', 'ycolor', [0.4 0.4 0.4], 'linewidth', 0.5, 'box', 'off');
axis(ax(2), 'square');
[~, pval2, ci, stats2] = ttest(y2(:, 1), y2(:, end));
bf10 = t1smpbf(stats2.tstat, 27);
fprintf('slope, bf10 = %f', bf10);
mysigstar(gca, [1 nbins], [0.9 0.9], pval, [0.4 0.4 0.4], 'down');

ax(1).YLabel.String = 'P(repeat)';
ax(2).YLabel.String = 'Slope';

switch whichMod
    case 'decision_pupil'
        xlabel('Previous trial pupil');
    case 'rt'
        xlabel('Previous trial RT');
    otherwise
        xlabel(sprintf('Previous trial %s', whichmodulator));
end
end
