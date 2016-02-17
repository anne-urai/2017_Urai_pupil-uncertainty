function fig5_psychFuncShift_Accuracy(whichmodulator, nbins, correctness)
global mypath;

if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('correctness', 'var'); correctness = []; end
if ~exist('nbins', 'var'); nbins = 3; end

warning('error', 'stats:glmfit:PerfectSeparation');

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
    case 'pupil-rt';
        whichMod = 'decision_pupil';
    case 'baselinepupil'
        whichMod = 'baseline_pupil';
    case 'evidence'
        whichMod = 'evidence';
end

subjects = 1:27;
clear grandavg;
lag = 1;

for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data(find(data.sessionnr > 1), :);
    
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
    
    % in this case, take out decision effects
    switch whichmodulator
        case 'fb-decpupil'
            data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        case 'pupil-rt',
            data.decision_pupil = projectout(data.decision_pupil, data.rt);
        case 'baselinepupil'
            data.baseline_pupil = circshift(data.baseline_pupil, -1);
        case 'evidence'
            data.evidence = abs(data.motionstrength);
    end
    
    % outcome vector need to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    
    % get an overall logistic fit for normalization
    grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, ...
        data.resp, 'binomial','link','logit');
    
    if isempty(correctness),
        % select for response and correctness
        trls = 1:size(data, 1);
    else
        trls = find(data.correct == correctness);
    end
    
    % fit an overall bias term
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
    
    %% split into pupil bins
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
        
        % take subset
        if ~isempty(correctness),
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
    
        % save betas
        grandavg.accuracy(sj, u, lag, :) = nanmean(thisdat.correct(thisdat.motionstrength < threshold));
        
    end % uncertainty bin
end % sj

% ========================================================= %
% plot
% ========================================================= %

subplot(441);
x = 1:nbins;

% Use HOLD and ERRORBAR, passing axes handles to the functions.
errorbar(x, mean(grandavg.accuracy), std(grandavg.accuracy) ./sqrt(27), ...
    '.-', 'color', 'k', 'markerfacecolor', 'k', 'markersize', 12);
%boundedline(x, mean(grandavg.accuracy), std(grandavg.accuracy) ./ sqrt(27), 'k');
axis tight; axis square;
xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high'});
ylabel('Next trial P(correct)');
box off; ylim([0.72 0.78]);
% low vs high
[~, pval] = ttest(grandavg.accuracy(:, 1), grandavg.accuracy(:, end));
mysigstar([1 nbins], [0.77 0.77], pval, 'k', 'down');

switch whichMod
    case 'decision_pupil'
        xlabel('Current trial pupil');
    case 'evidence'
        xlabel('Current trial evidence');
    case 'rt'
        xlabel('Current trial RT');
    otherwise
        xlabel(whichmodulator);
end

end
