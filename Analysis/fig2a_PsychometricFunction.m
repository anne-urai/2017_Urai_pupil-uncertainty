function fig2a_PsychometricFunction
% plot a logistic psychometric function for all subjects
% Anne Urai, 2015

global mypath;

subjects = 1:27;
% all the coherence levels that were used throughout
allcohs = [-0.3 -0.2 -0.1 -0.05 -0.025 -0.0125 -0.0063 0.0063 0.0125 0.025 0.05 0.1 0.2 0.3];
newx = linspace(min(allcohs), max(allcohs), 100);
cols = linspecer(5); cols = cols([1 4], :);

% preallocate
grandavg.betas = nan(length(subjects),  2);
grandavg.yfit = nan(length(subjects), 121);
grandavg.xpts = nan(length(subjects), length(allcohs));
grandavg.ypts = nan(length(subjects), length(allcohs));

for sj = unique(subjects),
    
    % disp(sj);
    % get all the data
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % make it such that we can plot the GA of datapoints
    LogisticFit = Psychometric_Logistic(data, 0, 0);
    for j = 1:length(allcohs),
        idx = find( abs(LogisticFit.coherences - allcohs(j)) < 0.0001);
        try
            grandavg.xpts(find(sj==subjects),  j) = LogisticFit.coherences(idx);
            grandavg.ypts(find(sj==subjects),  j) = LogisticFit.meanresp(idx);
        end
    end

    % do a quicker logistic fit without the lapse rate 
    x = data.coherence .* data.stim;
    y = data.resp; y(y==-1) = 0;
    
    % sort both
    [x, idx] = sort(x);
    y = y(idx);
    
    [beta,dev,stats] = glmfit(x, y, ...
        'binomial','link','logit');
    
    yfit = glmval(beta, newx,'logit');
    grandavg.betas(sj, :) = beta;
    
    % also save the curve itself
    grandavg.yfit(sj, :) = LogisticFit.y;
    
    
    if 0,
        % also get the individual d prime
        % use 0 (for weaker motion answer) instead of -1
        resp = data.resp; resp(resp==-1) = 0;
        stim = data.stim; stim(stim==-1) = 0;
        
        hit  = length(find(resp(find(stim == 1))==1)) / length(resp(find(stim == 1)));
        fa   = length(find(resp(find(stim == 0))==1)) / length(resp(find(stim == 0)));
        
        grandavg.dprime(sj) = norminv(hit) - norminv(fa);
        
        % see if there is a difference in the actual stimulus as well
        grandavg.difficulty(sj, :) = mean(abs(data.motionstrength));
    end
end

% plot the whole thing
plot(LogisticFit.x*100, nanmean(grandavg.yfit), 'color', [0.3 0.3 0.3]); % average function fit
hold on;

% then the datapoints, with standard deviation
p = ploterr(100*squeeze(nanmean(grandavg.xpts)), squeeze(nanmean(grandavg.ypts)), ...
    100*squeeze(nanstd(grandavg.xpts)), ...
    squeeze(nanstd(grandavg.ypts)) , 'k.', 'hhxy', 0.00000000000000000000000001);
p(1).MarkerSize = 8; p(1).Color = 'w';
box off; axis tight; set(gca, 'ytick', [0 0.5 1], 'xtick', [-30 -15 0 15 30]);
axis square;
xlabel('delta motion coherence'); ylabel('% response ''stronger''');
set(gca, 'xcolor', 'k', 'ycolor', 'k');

end


function [thisfit] = Psychometric_Logistic(dat, doBootstrap, grid)
%% fit a logistic psychometric function

% fprintf('this dataset has %d trials  ', height(dat));

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

% disp('fitting logistic...');

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

