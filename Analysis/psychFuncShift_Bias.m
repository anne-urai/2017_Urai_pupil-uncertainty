function b = psychFuncShift_Bias(whichmodulator, nbins, correctness, subjects)
global mypath;

if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('correctness', 'var'); correctness = []; end
if ~exist('nbins', 'var'); nbins = 3; end
lag = 1; % look at 1 trial in the past
if ~exist('subjects', 'var'); subjects = 1:27; end

warning('error', 'stats:glmfit:PerfectSeparation');

switch whichmodulator
    case 'pupil'
        whichMod = 'decision_pupil';
    case 'fbpupil'
        whichMod = 'feedback_pupil';
    case 'fb+decpupil'
        % see below for projecting out
        whichMod = 'feedback_pupil';
    case 'baseline_pupil'
        whichMod = 'baseline_pupil';
    otherwise
        whichMod = whichmodulator;
end

grandavg.logistic = nan(27, 2, nbins, 2);

% =========================================== %
% BIAS - ITERATE OVER TWO RESPONSES
% =========================================== %

for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data((data.sessionnr > 1), :); % get rid of residual learning effects
    
    switch whichmodulator
        case 'fb+decpupil'
            % in this case, take out decision effects
            data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        case 'pupil',
            data.decision_pupil = projectout(data.decision_pupil, zscore(log(data.rt + 0.1)));
        case 'baseline_pupil'
            % move to previous trial so that we split by current trial baseline
            data.baseline_pupil = circshift((data.baseline_pupil), -1);
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

% ========================================================= %
% plot
% ========================================================= %
%
x   = 1:nbins;
y1  = grandavg.logisticRep;

% Use HOLD and ERRORBAR, passing axes handles to the functions.
colors = cbrewer('qual', 'Set1', 8);
if isempty(correctness),
    thismarker = '.';
    thismarkersize = 14;
    thiscolor = [0 0 0];
else
    switch correctness
        case 1
            thismarker = '.';
            thismarkersize = 14;
            thiscolor = colors(2,:);
        case 0
            thismarker = '^';
            thismarkersize = 4;
            thiscolor = colors(1,:);
    end
end

% line to indicate 50 % repetition
plot([1 3], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;

h = ploterr(x, nanmean(y1), [], nanstd(y1) ./sqrt(27), ...
    'k-',  'abshhxy', 0);
set(h(1), 'color', thiscolor, ...
    'markersize', thismarkersize, ...
    'marker',thismarker);
if correctness == 0,
    set(h(1), 'markerfacecolor', 'w', 'markeredgecolor', thiscolor);
end
set(h(2), 'color', thiscolor); % line color

xticklabs = repmat({' '}, 1, nbins);
xticklabs{1} = 'low';
xticklabs{end} = 'high';
if nbins == 3, xticklabs{2} = 'med'; end

set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins, ...
    'xticklabel', xticklabs, 'ylim', [0.48 0.56], 'ytick', [0.5 0.56], ...
    'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'on');

if length(subjects) > 1,
    % low vs high
    [~, pval, ci, stats] = ttest(y1(:, 1), y1(:, end));
    bf10 = t1smpbf(stats.tstat, 27);
    fprintf('bias, t(%d) = %.3f, p = %.3f, bf10 = %f \n', stats.df, stats.tstat, pval, bf10);
    mysigstar(gca, [1 nbins], [0.55 0.55], pval, 'k', 'down');
else
    axis tight; xlim([0.5 nbins+0.5]);
end
axis square;
ylabel('P(repeat)');

switch whichMod
    case 'decision_pupil'
        xlabel('Previous trial pupil');
    case 'baseline_pupil'
        xlabel('Baseline pupil');
    case 'rt'
        xlabel('Previous trial RT');
    otherwise
        xlabel(sprintf('Previous trial %s', whichmodulator));
end

b = y1;
end
