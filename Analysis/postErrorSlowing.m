function [] = biasRT()
% check how pupil-linked uncertainty affects subsequent sensitivity, RT,
% side bias, absolute bias, or post-error slowing

global mypath;
nbins = 3;
clc;
clf;

% get all the data
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));

% discretize within each subject
data.prevPupilBins = nan(height(data), 1);
data.prev2RT       = nan(height(data), 1);

for sj = unique(data.subjnr)',
    data.prevPupilBins(sj == data.subjnr) = circshift(discretize( ...
        data.decision_pupil(sj == data.subjnr), [-inf quantile(data.decision_pupil(sj == data.subjnr), nbins-1) inf]), 1);
    
    % normalise RT
    data.prev2RT(sj == data.subjnr) = circshift(zscore(log(data.rt(sj == data.subjnr) + 0.1)), 2);
    data.rt(sj == data.subjnr)      = zscore(log(data.rt(sj == data.subjnr) + 0.1));
end

% dont include ones that cross
data.prevPupilBins([0; diff(data.trialnr)] ~= 1) = NaN;

% recode
data.resp(data.resp == -1) = 0;

% apply glm to get overall bias for each bin
myFun = @(x, y) transpose(glmfit(x, [y ones(size(y))], 'binomial','link','logit'));
grandavg = rowfun(myFun, data, 'inputvariables', {'motionstrength', 'resp'}, ...
    'groupingvariables', {'subjnr', 'prevPupilBins'}, 'outputvariablenames', 'beta');
grandavg.bias       = grandavg.beta(:, 1);
grandavg.absbias    = abs(grandavg.beta(:, 1));
grandavg.slope      = grandavg.beta(:, 2);

% next trial RT
grandavg = join(grandavg, rowfun(@nanmean, data, 'inputvariables', {'rt'}, ...
    'groupingvariables', {'subjnr', 'prevPupilBins'}, 'outputvariablenames', {'RT'}));

% next trial accuracy
grandavg = join(grandavg, rowfun(@nanmean, data, 'inputvariables', {'correct'}, ...
    'groupingvariables', {'subjnr', 'prevPupilBins'}, 'outputvariablenames', {'accuracy'}));

% post-error slowing, see Dutilh et al. 2012
data.rt = data.rt - data.prev2RT; % for this trial's RT, normalize by the RT 2 trials in the past

% get data where the last trial was an error, the trial before that was
% correct, and the current trial is correct again. That means we're
% subtracting the 2nd last trial's RT, with only the error in between, and
% taking the mean of the current post-error RT.
errtrls = (data.correct == 1 & circshift(data.correct, 1) == 0 & circshift(data.correct, 2) == 1);
grandavg = join(grandavg, rowfun(@nanmean, data(errtrls, :), 'inputvariables', {'rt'}, ...
    'groupingvariables', {'subjnr', 'prevPupilBins'}, 'outputvariablenames', {'slowingError'}), ...
    'keys', {'subjnr', 'prevPupilBins'});

corrtrls = (data.correct == 1 & circshift(data.correct, 1) == 1 & circshift(data.correct, 2) == 1);
grandavg = join(grandavg, rowfun(@nanmean, data(corrtrls, :), 'inputvariables', {'rt'}, ...
    'groupingvariables', {'subjnr', 'prevPupilBins'}, 'outputvariablenames', {'slowingCorrect'}), ...
    'keys', {'subjnr', 'prevPupilBins'});

%% plot all of those
plotVars = {'bias', 'absbias', 'RT', 'accuracy', 'slowingError'};
for p = 1:length(plotVars),
    subplot(3,3,p);
    errorbar(1:nbins, splitapply(@mean, grandavg.(plotVars{p}), findgroups(grandavg.prevPupilBins)), ...
        splitapply(@sem, grandavg.(plotVars{p}), findgroups(grandavg.prevPupilBins)), '.k-', 'markersize', 12);
    
    % ttest between lowest and highest bins
    [~, pval, ~, stats] = ttest(grandavg.(plotVars{p})((grandavg.prevPupilBins == 1)), ....
        grandavg.(plotVars{p})(grandavg.prevPupilBins == nbins));
    
    xlim([0.5 nbins+0.5]); box off;
    mysigstar(gca, [1 nbins], max(get(gca, 'ylim')), pval);
    xlabel('Previous pupil bin');
    ylabel(plotVars{p});
    
    % check for no difference
    bf10 = t1smpbf(stats.tstat, 27);
    assert(bf10 < 0.3, '%s changes with pupil!', plotVars{p});
end
disp('nothing changes with pupil');

end
