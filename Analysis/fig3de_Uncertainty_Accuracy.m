function [grandavg] = fig3de_Uncertainty_Accuracy(nbins)
% plots uncertainty by accuracy both for the modelfits and the pupil

if ~exist('nbins', 'var'), nbins = 100; end
close all;

subjects      = 1:27;
fitIndividual = false;

for sj = unique(subjects),
    
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    % divide into bins
    [ grandavg.pup(sj, :), grandavg.acc(sj, :), stdx, stdy] = ...
        divideintobins(data.decision_pupil, data.correct, nbins);
    
    % do a logistic regression
    x = [zscore(abs(data.motionstrength)) zscore(data.rt) zscore(data.decision_pupil) ];
    
    x(data.trialnr < 2, :) = nan;
    
    [b,dev,stats] = glmfit(x, data.correct, ...
        'binomial','link','logit');
    % save for later
    grandavg.b(sj, :) = b;
    
    if fitIndividual,
        subplot(5,6,sj);
        yfit = glmval(b, x,'logit');
        hold on;
        plot(grandavg.pup(sj, :), grandavg.acc(sj,:), 'b.');
        axis tight; ylim([0.5 1]); xlims = get(gca, 'xlim');
        plot(x,yfit,'k');
        xlim(xlims);
        title(sprintf('P%02d', sj));
    end
    
end

% now, plot
cols = linspecer(3);
figure;
subplot(441);

grandavg.acc = 100 * grandavg.acc;
boundedline(1:nbins, nanmean(grandavg.acc), ...
    nanstd(grandavg.acc) / sqrt(length(subjects)), 'cmap', [0 0 0]);

xlabel('Pupil response');
ylabel({'Percent correct'});
axis tight; axis square;
set(gca, 'ytick', [65 75 85]);
set(gca, 'xtick', [1  nbins/2 nbins], 'xticklabel', {'low', 'medium', 'high'});
% offsetAxes;
box off; ylim([65 85]);

if 0,
    subplot(463);
    
    [~, pval(1)] = permtest(grandavg.betas(:, 1));
    [~, pval(2)] = permtest(grandavg.betas(:, 2));
    [~, pval(3)] = permtest(grandavg.betas(:, 3));
    [~, pval(4)] = permtest(grandavg.betas(:, 4));
    
    h = ploterr(1:size(grandavg.b, 2), squeeze(mean(grandavg.b)), [], squeeze(std(grandavg.b)) / sqrt(length(subjects)), ...
        'k.', 'hhxy', 0.0000000000000001);
    set(h(1), 'marker', 'none');
    hold on;
    bar(1:size(grandavg.b, 2), mean(grandavg.b), 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none', 'barwidth', 0.4);
    sigstar({[1 1], [2 2], [3 3], [4 4]}, pval);
    
    set(gca, 'xtick', 1:4, 'xticklabel', {'Intercept', '|\Delta|', 'RT', 'Pupil'}, 'xticklabelrotation', -40)
    box off;
    ylabel('Beta weights (a.u.)');
    set(gca, 'xcolor', 'k', 'ycolor', 'k')
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig3d_pupil_accuracy.pdf'));

%% fit logistic slope on high vs low pupil bins
subjects = 1:27;
allcohs = [-0.3 -0.2 -0.1 -0.05 -0.025 -0.0125 -0.0063 0.0063 0.0125 0.025 0.05 0.1 0.2 0.3];

grandavg.xpts  = nan(length(subjects), 2, length(allcohs));
grandavg.ypts  = nan(length(subjects), 2, length(allcohs));
grandavg.betas = nan(length(subjects), 2, 2);

for sj = unique(subjects),
    
    disp(sj);
    % get all the data
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    % divide into low and high pupil bins
    puptrls{1} = find(data.decision_pupil < quantile(data.decision_pupil, 0.5));
    puptrls{2} = find(data.decision_pupil > quantile(data.decision_pupil, 0.5));
    
    % clear puptrls; % puptrls{1} = 1:height(data);
    
    for b = 1:length(puptrls),
        thisdat = data(puptrls{b}, :);
        
        % subplot(5,6,sj);
        LogisticFit = Psychometric_Logistic(thisdat, 0, 0);
        
        % make it such that we can plot the GA
        for j = 1:length(allcohs),
            idx = find( abs(LogisticFit.coherences - allcohs(j)) < 0.0001);
            try
                grandavg.xpts(find(sj==subjects), b, j) = LogisticFit.coherences(idx);
                grandavg.ypts(find(sj==subjects), b, j) = LogisticFit.meanresp(idx);
            end
        end
        
        % fit on all the trials, dont average beforehand!
        x = thisdat.motionstrength;
        y = thisdat.resp; y(y==-1) = 0;
        
        % sort both
        [x, idx] = sort(x);
        y = y(idx);
        
        [beta,dev,stats] = glmfit(x, [y ones(size(y))], ...
            'binomial','link','logit');
        newx = -3:0.01:3;
        yfit = glmval(beta, newx,'logit');
        % plot(newx, yfit,'-', 'color', cols(b, :));
        
        % also save the curve itself
        grandavg.betas(sj, b, :) = beta;
        grandavg.yfit(sj, b, :) = yfit;
        
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
subplot(4,5,3);
h = ploterr(1:2, squeeze(nanmean(grandavg.betas(:, :, 2))), [], ...
    squeeze(std(grandavg.betas(:, :, 2))) / sqrt(length(subjects)), 'k.', 'hhxy', 0.00001);
set(h(1), 'marker', 'none');

hold on;
bar(1, mean(grandavg.betas(:, 1, 2)), 'facecolor', [0.6 0.6 0.6], 'edgecolor', 'none', 'barwidth', 0.4);
bar(2, mean(grandavg.betas(:, 2, 2)), 'facecolor', [0.4 0.4 0.4], 'edgecolor', 'none', 'barwidth', 0.4);
% do ttest on regression coefficients
[~, pval(3), ~, stat] = ttest(grandavg.betas(:, 1, 1), grandavg.betas(:, 1, 2));
ylim([0.9 1.5]); xlim([0.5 2.5]); box off;
s1 = sigstar({[1 2]}, pval(3), 0); ylim([0.9 1.5]); set(gca, 'ytick', [1 1.5]);
xlabel('Pupil response'); set(gca, 'xtick', 1:2, 'xticklabel', {'low', 'high'});
ylabel('Perceptual sensitivity');

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig3e_pupil_sensitivity.pdf'));
savefast('~/Data/pupilUncertainty/GrandAverage/grandavg_logistic_bypupil.mat', 'grandavg', 'allcohs', 'subjects', 'cols');

%%
% !!! rather than doing a permutation test on the logistic slope
% coefficients, run a proper mixed effects logistic regression
% fixed effects: slope
if 0,
    clear; clc;
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_allsj.csv'));
    
    % normalize the pupil response for each subject
    % this way, we can interpret the coefficient of the pupil fixed effect in
    % units of standard deviation
    for sj = unique(data.subjnr)',
        data.decision_pupil(sj==data.subjnr) = zscore(data.decision_pupil(sj==data.subjnr));
    end
    
    data.absmotion = abs(data.motionstrength);
    
    mdl = fitglme(data, 'correct ~ 1 + decision_pupil + absmotion + (1|subjnr) + (1|sessionnr) ', ...
        'Distribution', 'Binomial', 'Link', 'Logit');
    disp(mdl);
end

%% compute ROC AUC based on pupil
clear all; clc; close all;

subjects = 1:27;
for sj = unique(subjects),
    
    % get all the data
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    out = rocAnalysis(data.decision_pupil(data.correct==1), ...
        data.decision_pupil(data.correct==0), 0, 1);
    
    grandavg.roc(sj)    = out.i;
    grandavg.pval(sj)   = out.p;
end
savefast('~/Data/pupilUncertainty/GrandAverage/pupilCorrectnessROC.mat', 'grandavg');

load('~/Data/pupilUncertainty/GrandAverage/pupilCorrectnessROC.mat');
% test across the group, roc auc values are normally distributed
[h, pval ci, stats] = ttest(grandavg.roc, 0.5, 'tail', 'both')

end

function [thisfit] = Psychometric_Logistic(dat, doBootstrap, grid)
%% fit a logistic psychometric function

fprintf('this dataset has %d trials  ', height(dat));

% define the functions used for fitting
logistic = @(coh, bias, sensitivity, lower, upper) lower+(1-lower-upper)*(1./(1+exp(-1*(sensitivity).*(coh-bias))));
% as in Kahnt et al, but with additional lower and upper lapse rates

% loglikelihood of this
LL_logistic = @(intensities, results, bias, sensitivity, lower, upper) -sum(results.*log(logistic(intensities, ...
    bias, sensitivity, lower, upper)) + (1-results).*log(1-logistic(intensities, bias, sensitivity, lower, upper)));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE DATAPOINTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add a col to make my life easier
stim = dat.stim .* dat.coherence;
resp = dat.resp; resp(resp == -1) = 0;

coherences = unique(stim)';
for icoh = 1:length(coherences),
    meanresp(icoh) = mean(resp(find(dat.stim==sign(coherences(icoh)) & ...
        dat.coherence == abs(coherences(icoh)))));
end

%assert(~any(isnan(meanresp)), 'nans in mean response curves');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT LOGISTIC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('fitting logistic...');

% function evaluation params
options                 = optimset('Display', 'notify') ;
options.MaxFunEvals     = 100000000000000000000000000;
options.MaxIter         = 500000;
options.TolX            = 0.00000001;
options.TolFun          = 0.00000001;
options.Robust          = 'on';

% initial values for fminsearch
guessbeta(1) = 0; % bias
guessbeta(2) = 15; % sensitivity, slope
guessbeta(3) = 0; % lower rate
guessbeta(4) = 0; % upper rate (lapse rate)

% lower bound
lowerbound(1)   = -.2;
lowerbound(2)   = 0;
lowerbound(3)   = 0;
lowerbound(4)   = 0;

% upper bound
upperbound(1)   = .2;
upperbound(2)   = 100; % what is the maximum value this could take? be liberal!
upperbound(3)   = .1; % 10% lapse seems sufficient
upperbound(4)   = .1;

if grid,
    
    cnt = 1;
    severalbetas = nan((numel(lowerbound(1):0.1:upperbound(1))+...
        numel(lowerbound(2):upperbound(2))+numel(lowerbound(4):0.1:upperbound(4))), 4);
    severalfits_fval = nan((numel(lowerbound(1):0.1:upperbound(1))+...
        numel(lowerbound(2):upperbound(2))+numel(lowerbound(4):0.1:upperbound(4))), 1);
    
    % LOOP THROUGH A GRID OF STARTING POINTS
    for bias = lowerbound(1):0.1:upperbound(1),
        for slope = lowerbound(2):upperbound(2),
            for lapse = lowerbound(3):0.1:upperbound(3),
                
                guessbeta(1) = bias;
                guessbeta(2) = slope;
                guessbeta(3) = lapse;
                guessbeta(4) = lapse;
                
                % find optimal values for beta using fminsearch
                [severalbetas(cnt, :), severalfits_fval(cnt)] = fminsearchbnd(@(beta) LL_logistic(coherences, ...
                    meanresp, beta(1), beta(2), beta(3), beta(4)), guessbeta, lowerbound, upperbound, options);
                cnt = cnt + 1;
            end
        end
        
    end
    
    % from all the tried starting points, select the best fit
    [val, bestfit]      = min(severalfits_fval);
    thisfit.mean.beta   = severalbetas(bestfit, :);
    thisfit.mean.fval   = val;
    fprintf('betas: %f %f %f %f, fval: %f', thisfit.mean.beta, thisfit.mean.fval);
    
else
    
    [thisfit.mean.beta, thisfit.mean.fval] = fminsearchbnd(@(beta) LL_logistic(stim, ...
        resp, beta(1), beta(2), beta(3), beta(4)), [0 15 0 0], lowerbound, upperbound, options);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NONPARAMETRIC BOOTSTRAP FOR THE DATAPOINTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doBootstrap,
    
    % BOOTSTRAP PARAMS
    nboot               = 1000;
    boot_resp           = nan(length(coherences), nboot);
    
    disp('nonparametric bootstrapping...');
    tic;
    x = linspace(min(coherences), max(coherences));
    
    for iboot = 1:nboot,
        for icoh = 1:length(coherences),
            
            % get a random sample from all the datapoints with replacement - only for this level of difficulty
            sampleidx   = datasample(find(stim==coherences(icoh)), ...
                length(find(stim==coherences(icoh))), 'Replace', true);
            boot_resp(icoh, iboot)    = mean(resp(sampleidx));
        end
        
        % FIT TO THESE BOOTSTRAPPED DATA, TO GET A CONFIDENCE INTERVAL AROUND THE SLOPE
        [thisbeta(iboot, :), ~] = fminsearchbnd(@(beta) LL_logistic(coherences, ...
            boot_resp(:, iboot)', beta(1), beta(2), beta(3), beta(4)), thisfit.mean.beta, lowerbound, upperbound, options);
    end
    toc;
    
    % compute confidence intervals as percentiles of the distribution
    CI.lower        = prctile(boot_resp', 2.5);
    CI.upper        = prctile(boot_resp', 97.5);
    
    % also get confidence intervals around the parameters
    CI.betas.lower = prctile(thisbeta, 2.5);
    CI.betas.upper = prctile(thisbeta, 97.5);
    
else
    CI.lower = nan(size(coherences));
    CI.upper = nan(size(coherences));
    CI.betas.lower = nan(size(thisfit.mean.beta));
    CI.betas.upper = nan(size(thisfit.mean.beta));
    
end

% get the x and y values for plotting
thisfit.x = -.3:0.005:0.3;
thisfit.y = logistic(thisfit.x, thisfit.mean.beta(1), thisfit.mean.beta(2), thisfit.mean.beta(3), thisfit.mean.beta(4));

thisfit.CI = CI;
thisfit.meanresp = meanresp;
thisfit.coherences = coherences;

end
