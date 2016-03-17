function psychFuncShift_Bias_Slope(whichmodulator, nbins, correctness)
global mypath;

if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('correctness', 'var'); correctness = []; end
if ~exist('nbins', 'var'); nbins = 3; end

warning('error', 'stats:glmfit:PerfectSeparation');

switch whichmodulator
    case 'pupil'
        whichMod = 'decision_pupil';
    case 'fbpupil'
        whichMod = 'feedback_pupil';
    case 'fb+decpupil'
        % see below for projecting out
        whichMod = 'feedback_pupil';
    case 'rt'
        whichMod = 'rt';
end


% =========================================== %
% BIAS - ITERATE OVER TWO RESPONSES
% =========================================== %

subjects = 1:27;
clear grandavg;
lag = 1;

for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data(find(data.sessionnr > 1), :);
    
    % in this case, take out decision effects
    switch whichmodulator
        case 'fb+decpupil'
            data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        case 'pupil',
            data.decision_pupil = projectout(data.decision_pupil, zscore(log(data.rt + 0.1)));
        case 'rt'
            data.rt = projectout(zscore(log(data.rt+0.1)), data.decision_pupil);
    end
    
    % outcome vector need to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    
    % get an overall logistic fit for normalization
    grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, ...
        data.resp, 'binomial','link','logit');
    
    % previous response
    resps = [0 1];
    for r = 1:2,
        
        % split into pupil bins
        if nbins > 2,
            if ~isempty(correctness),
                uncQs = quantile(data.(whichMod)(data.resp == resps(r) & data.correct == correctness), nbins - 1);
            else
                uncQs = quantile(data.(whichMod)(data.resp == resps(r)), nbins - 1);
            end
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
            
            try
                b = glmfit(thisdat.motionstrength, thisdat.resp, ...
                    'binomial','link','logit');
            catch
                warning('putting nan in betas!');
                b = nan(size(b));
            end
            
            % save betas
            grandavg.logistic(sj, r, u, lag, :) = b;
            
        end % uncertainty bin
    end % resp
end % sj

% ========================================================= %
% normalize by their overall bias
% ========================================================= %

grandavg.logistic = bsxfun(@minus, grandavg.logistic, grandavg.overallLogistic(:, 1));

% ========================================================= %
% transform to probability saying '1'
% ========================================================= %

grandavg.logistic = exp(grandavg.logistic)./(1+exp(grandavg.logistic));

% ========================================================= %
% make one plot with repetition (rather than response) bias
% ========================================================= %

resp1 = -1 * (squeeze(grandavg.logistic(:, 1, :, :, 1)) - 0.5) + 0.5;
resp2 = squeeze(grandavg.logistic(:, 2, :, :, 1));

% since this is centred at 0.5, treat it that way
grandavg.logisticRep = (resp1 + resp2) ./ 2;

% ========================================================= %
% combine the lag groups
% ========================================================= %

grandavg.logisticRep = squeeze(nanmean(grandavg.logisticRep(:, :, 1), 3));
grandavg.logistic = squeeze(nanmean(grandavg.logistic(:, :, :, 1, :), 4));
grandavgBias = grandavg;

% =========================================== %
% SLOPE - DO NOT ITERATE OVER TWO RESPONSES   %
% =========================================== %

for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data(find(data.sessionnr > 1), :);
    
    % in this case, take out decision effects
    % in this case, take out decision effects
    switch whichmodulator
        case 'fb+decpupil'
            data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        case 'pupil',
            data.decision_pupil = projectout(data.decision_pupil, zscore(log(data.rt + 0.1)));
        case 'rt'
            data.rt = projectout(zscore(log(data.rt+0.1)), data.decision_pupil);
    end
    
    % outcome vector needs to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    
    % previous response
    if nbins > 2,
        if ~isempty(correctness),
            uncQs = quantile(data.(whichMod)(data.correct == correctness), nbins - 1);
        else
            uncQs = quantile(data.(whichMod), nbins - 1);
        end
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
        trls(removeTrls) = [];
        
        % fit logistic regression
        thisdat = data(laggedtrls, :);
        grandavg.logistic(sj, u, lag, :) = glmfit(thisdat.motionstrength, thisdat.resp, ...
            'binomial','link','logit');
        
    end % uncertainty bin
end % sj

% replace
grandavg.logisticSlope = squeeze(mean(grandavg.logistic(:, :, 1, 2), 3));
grandavgSlope = grandavg;

% ========================================================= %
% plot
% ========================================================= %

x = 1:nbins;
y1 = grandavgBias.logisticRep;
y2 = grandavgSlope.logisticSlope;

% Store the axes handles produced by PLOTYY
[ax, p1, p2] = plotyy(x, mean(y1), x, mean(y2));
p1.Color = 'k';
p2.Color = [0.4 0.4 0.4];

% Use HOLD and ERRORBAR, passing axes handles to the functions.
hold(ax(1), 'on');
errorbar(ax(1), x, mean(y1), std(y1) ./sqrt(27), ...
    '.-', 'color', 'k', 'markerfacecolor', 'k', 'markersize', 12);

set(ax(1), 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins, ...
    'xticklabel', {'low', 'med', 'high'}, 'ylim', [0.48 0.56], 'ytick', [0.49 0.51 0.53], ...
    'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off');
% low vs high
[~, pval, ci, stats] = ttest(y1(:, 1), y1(:, end));
disp(stats);
disp(pval);
mysigstar([1 nbins], [0.53 0.53], pval, 'k', 'down');
axis(ax(1), 'square');

% second errorbar, grey
hold(ax(2), 'on');
errorbar(ax(2), x, mean(y2), std(y2) ./sqrt(27), ...
    '.-', 'color', [0.4 0.4 0.4], 'markerfacecolor', [0.4 0.4 0.4], 'markersize', 12);

set(ax(2), 'xlim', [0.5 nbins+0.5], 'ylim', [0.40 0.9], 'ytick', [0.8 0.9], ...
    'xcolor', 'k', 'ycolor', [0.4 0.4 0.4], 'linewidth', 0.5, 'box', 'off');
axis(ax(2), 'square');
[~, pval] = ttest(y2(:, 1), y2(:, end));
mysigstar([1 nbins], [0.9 0.9], pval, [0.4 0.4 0.4], 'down');

ax(1).YLabel.String = 'P(repeat)';
ax(2).YLabel.String = 'Slope';

switch whichMod
    case 'decision_pupil'
        xlabel('Previous trial pupil');
    case 'rt'
        xlabel('Previous trial RT');
    otherwise
        xlabel(whichmodulator);
end

end
