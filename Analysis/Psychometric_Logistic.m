
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
