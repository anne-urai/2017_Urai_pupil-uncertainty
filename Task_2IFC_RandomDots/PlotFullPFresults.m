function [datapoints1, datapoints2, datapoints3, fit1, fit2, fit3] = PlotFullPFresults(dataset, fitting, bootstrapping, BAC)%% fit a psychometric function

load(dataset);
warning off;

if ~exist('fitting')
    fitting = 'fit';
end

if ~exist('bootstrapping')
    bootstrapping = 'bootstrap';
end

% define cumulative Weibull including chance level of 0.5
cumWB = @(coh, threshold, slope, lower, upper) lower+(1- lower -upper)*(1-exp(-1*(coh./threshold).^slope));

%define logLikelihood of the cumulative WB
LL_cumWB = @(intensities, results, threshold, slope, lower, upper) -sum(results.*log(cumWB(intensities, ...
    threshold, slope, lower, upper)) + (1-results).*log(1-cumWB(intensities, threshold, slope, lower, upper)));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSYCHOMETRIC FUNCTION NR1 - RAW
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reshape
datapoints1.all.difficulty  = reshape(setup.difficulty', 1, setup.nblocks*setup.ntrials);
datapoints1.all.correct     = reshape(results.correct', 1, setup.nblocks*setup.ntrials);

% average over difficulty levels
datapoints1.mean.difficulty = unique(datapoints1.all.difficulty);
for idiff = 1:length(datapoints1.mean.difficulty),
    datapoints1.mean.correct(idiff) = nanmean(datapoints1.all.correct(datapoints1.all.difficulty==datapoints1.mean.difficulty(idiff)));
end

% ------- BOOTSTRAP DATAPOINTS ---------------------------------------------------------------------------

switch bootstrapping,
    case 'bootstrap',
        
        % BOOTSTRAP PARAMS
        nboot = 500;
        
        % preallocate
        datapoints1.bootstrap.correct = nan(length(datapoints1.mean.difficulty), nboot);
        
        disp('bootstrapping...');
        
        % NONPARAMETRIC bootstrap the obtained values
        for iboot = 1:nboot,
            % do this separately per level of difficulty!
            for idiff = 1:length(datapoints1.mean.difficulty),
                % get a random sample from all the datapoints with replacement - only for this level of difficulty
                [bootstrappeddata] = datasample(datapoints1.all.correct(datapoints1.all.difficulty==datapoints1.mean.difficulty(idiff)), ...
                    length(datapoints1.all.difficulty(datapoints1.all.difficulty==datapoints1.mean.difficulty(idiff))), 'Replace', true);
                
                datapoints1.bootstrap.correct(idiff, iboot) = mean(bootstrappeddata);
            end
        end
        
    case 'nobootstrap',
        disp('not bootstrapping CIs');
end

% ------- FIT CUMULATIVE WEIBULL ---------------------------------------------------------------------------

switch fitting,
    case 'fit',
        
        % function evaluation params
        options.MaxFunEvals = 5000000;
        options.MaxIter = 500000;
        options.TolX = 0.00000001;
        options.TolFun = 0.00000001;
        options = optimset('Display', 'off') ;
        options.Robust = 'on';
        
        % initial values for fminsearch
        guessbeta(1) = .15; %threshold
        guessbeta(2) = 2; %slope
        guessbeta(3) = 0.5; %lower rate
        guessbeta(4) = .01; % upper rate (lapse rate)
        
        % lower bound of params
        lowerbound(1)   = 0;
        lowerbound(2)   = 0;
        lowerbound(3)   = .5;
        lowerbound(4)   = 0;
        
        % upper bound of params
        upperbound(1)   = 1;
        upperbound(2)   = +inf;
        upperbound(3)   = .5;
        upperbound(4)   = .5;
        
        fit1.CI = [2.5 97.5]; %95% confidence intervals
        
        % find optimal values for beta using fminsearch
        if exist('simulannealbnd.m') == 2,
            % check if the global optimization toolbox is present
            [fit1.mean.beta, fit1.mean.fval] = simulannealbnd(@(beta) LL_cumWB(datapoints1.mean.difficulty, ...
                datapoints1.mean.correct, beta(1), beta(2), beta(3), beta(4)), guessbeta, lowerbound, upperbound, options);
            % simulated annealing finds a global minimum so is more
            % accurate than fminsearch
        else
            [fit1.mean.beta, fit1.mean.fval] = fminsearchbnd(@(beta) LL_cumWB(datapoints1.mean.difficulty, ...
                datapoints1.mean.correct, beta(1), beta(2), beta(3), beta(4)), guessbeta, lowerbound, upperbound, options);
        end
        
        % ----------- PARAMETRIC BOOTSTRAP FOR WEIBULL PARAMETER CIs -------
        
        switch bootstrapping
            case 'bootstrap',
                for iboot = 1:nboot,
                    
                    % get the percent correct based on the weibull fit
                    prob = cumWB(datapoints1.all.difficulty, fit1.mean.beta(1), fit1.mean.beta(2), fit1.mean.beta(3), fit1.mean.beta(4));
                    % create a mock response set
                    % this is the source of randomness in the bootstrap
                    boot_resp = binornd(ones(1,size(datapoints1.all.correct,2)), prob(1:size(datapoints1.all.correct,2)));
                    
                    % change this into parameters for weibull fit
                    for idiff = 1:length(datapoints1.mean.difficulty),
                        thisdiff = datapoints1.mean.difficulty(idiff);
                        % parametric bootstrap for weibull parameters
                        boot_corr(idiff) = nanmean(length(find(datapoints1.all.difficulty == thisdiff & boot_resp == 1))/length(find(datapoints1.all.difficulty==thisdiff)));
                        
                    end % coh
                    
                    % fit the weibull function to this set of betas
                    [fit1.boot.beta(:,iboot), fit1.boot.fval(iboot)] = fminsearchbnd(@(beta) LL_cumWB(datapoints1.mean.difficulty, ...
                        boot_corr, beta(1), beta(2), beta(3), beta(4)), guessbeta, lowerbound, upperbound, options);
                    
                end % boot
                
            case 'nobootstrap'
                %  do nothing
        end
end

% ------- PLOT ---------------------------------------------------------------------------

figure; %set(gcf,'units','normalized','outerposition', [0 0 1 1], 'color', 'w');% maximize the figure
set(gcf,'PaperOrientation','landscape');

subplot(221);

switch fitting,
    case 'fit',
        % get the x and y values for plotting
        fit1.x = linspace(min(datapoints1.mean.difficulty),max(datapoints1.mean.difficulty));
        % fit the cumulative Weibull
        fit1.y = cumWB(fit1.x, fit1.mean.beta(1), fit1.mean.beta(2), fit1.mean.beta(3), fit1.mean.beta(4));
        
        %plot
        semilogx(fit1.x*100, fit1.y, '-k', 'LineWidth', 2); hold on;
        % add weibull plots for the confidence intervals
        
        % add param info
        text(20, .5, sprintf('threshold %.2f %% \n slope %.2f \n guessrate %.2f %% \n lapse rate %.2f %%', ...
            fit1.mean.beta(1)*100, fit1.mean.beta(2), fit1.mean.beta(3)*100, fit1.mean.beta(4)*100));
end

% plot the datapoints by themselves
semilogx(datapoints1.mean.difficulty*100, datapoints1.mean.correct, 'bo', 'MarkerSize', 8,  'MarkerFaceColor', 'b'); hold on;
hold on;
% add a thin dotted line to indicate alpha
fit1.correct = fit1.y(dsearchn(fit1.x', fit1.mean.beta(1)));

switch bootstrapping,
    case 'bootstrap',
        
        %  the confidence interval errorbars for the datapoints
        datapoints1.bootstrap.CI.lower = prctile(datapoints1.bootstrap.correct', fit1.CI(1));
        datapoints1.bootstrap.CI.upper = prctile(datapoints1.bootstrap.correct', fit1.CI(2));
        errorbar(datapoints1.mean.difficulty*100, datapoints1.mean.correct, ...
            datapoints1.mean.correct - datapoints1.bootstrap.CI.lower, datapoints1.bootstrap.CI.upper - datapoints1.mean.correct, 'b.');
        
        % confidence intervals for the beta values
        fit1.boot.CI.lower = prctile(fit1.boot.beta(1,:)', fit1.CI(1));
        fit1.boot.CI.upper = prctile(fit1.boot.beta(1,:)', fit1.CI(1));
        
end

title(sprintf('Raw performance'));
xlabel('Motion coherence from baseline (%)'); ylabel('Accuracy (% correct)');
ylim([.39 1.01]);
set(gca, 'XTick', datapoints1.mean.difficulty*100);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSYCHOMETRIC FUNCTION NR2 - D' CORRECTED
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(222);

% reshape
datapoints2.all.coherence       = reshape(dots.coherence', 1, setup.nblocks*setup.ntrials);
datapoints2.all.increment       = reshape(setup.increment', 1, setup.nblocks*setup.ntrials);
datapoints2.all.response        = reshape(results.response', 1, setup.nblocks*setup.ntrials);

% use 0 (for weaker motion answer) instead of -1
datapoints2.all.response(datapoints2.all.response == -1) = 0;
%datapoints2.all.increment(datapoints2.all.increment == 270) = 0;

% ------- BOOTSTRAP DATAPOINTS ---------------------------------------------------------------------------

switch bootstrapping,
    case 'bootstrap',
        % BOOTSTRAP PARAMS
        nboot = 500;
        
        disp('bootstrapping...');
        % bootstrap the obtained values
        for iboot = 1:nboot,
            % do this separately per level of difficulty!
            
            for idiff = 1:length(datapoints1.mean.difficulty),
                
                % find the index of the samples with this level of
                % difficulty
                thisdiffidx = find(datapoints1.all.difficulty==datapoints1.mean.difficulty(idiff));
                % get a random sample from all the datapoints with replacement - only for this level of difficulty
                sampleidx = datasample(thisdiffidx, length(thisdiffidx), 'Replace', true);
                
                % get the data for this bootstrap
                datapoints2.this.response   = datapoints2.all.response(sampleidx);
                datapoints2.this.difficulty = datapoints1.all.difficulty(sampleidx);
                datapoints2.this.increment  = datapoints2.all.increment(sampleidx);
                
                % for each difficulty level, calculate the hit and FA rate
                datapoints2.bootstrap.hit(idiff, iboot)  = length(find(datapoints2.this.response(find(datapoints2.this.difficulty == datapoints1.mean.difficulty(idiff) & datapoints2.this.increment == 1))==1)) / length(datapoints2.this.response(find(datapoints2.this.difficulty == datapoints1.mean.difficulty(idiff) & datapoints2.this.increment == 1)));
                datapoints2.bootstrap.fa(idiff, iboot)   = length(find(datapoints2.this.response(find(datapoints2.this.difficulty == datapoints1.mean.difficulty(idiff) & datapoints2.this.increment == -1))==1)) / length(datapoints2.this.response(find(datapoints2.this.difficulty == datapoints1.mean.difficulty(idiff) & datapoints2.this.increment == -1)));
                
                % from email Tobi 10 Feb
                datapoints2.bootstrap.dprime(idiff,iboot)     = norminv(datapoints2.bootstrap.hit(idiff, iboot), 0, 1) - norminv(datapoints2.bootstrap.fa(idiff, iboot), 0, 1);
                datapoints2.bootstrap.correctedy(idiff,iboot) = normcdf(0.5*(datapoints2.bootstrap.dprime(idiff,iboot)));
                
            end
            
        end
    case 'nobootstrap',
        disp('not bootstrapping CIs');
end

% average over the bootstrapped values for the mean and CIs
switch bootstrapping,
    case 'bootstrap'
        datapoints2.mean.dprime         = mean(datapoints2.bootstrap.dprime, 2);
        datapoints2.mean.correctedy     = mean(datapoints2.bootstrap.correctedy, 2)';
        
    case 'nobootstrap'
        for idiff = 1:length(datapoints1.mean.difficulty),
            
            % for each difficulty level, calculate the hit and FA rate
            datapoints2.mean.hit(idiff)  = length(find(datapoints2.all.response(find(datapoints1.all.difficulty == datapoints1.mean.difficulty(idiff) & datapoints2.all.increment == 1))==1)) / length(datapoints2.all.response(find(datapoints1.all.difficulty == datapoints1.mean.difficulty(idiff) & datapoints2.all.increment == 1)));
            datapoints2.mean.fa(idiff)   = length(find(datapoints2.all.response(find(datapoints1.all.difficulty == datapoints1.mean.difficulty(idiff) & datapoints2.all.increment == 0))==1)) / length(datapoints2.all.response(find(datapoints1.all.difficulty == datapoints1.mean.difficulty(idiff) & datapoints2.all.increment == 0)));
        end
        
        % from email Tobi 10 Feb
        datapoints2.mean.dprime              = norminv(datapoints2.mean.hit) - norminv(datapoints2.mean.fa);
        datapoints2.mean.correctedy          = normcdf(0.5*(datapoints2.mean.dprime));
end

% ------- FIT A CUMULATIVE WEIBULL ---------------
switch fitting,
    case 'fit',
        
        % !!! use the values that were already stored for the previous function fit
        
        % find optimal values for beta using fminsearch
        [fit2.mean.beta, fit2.mean.fval] = fminsearchbnd(@(beta) LL_cumWB(datapoints1.mean.difficulty, ...
            datapoints2.mean.correctedy, beta(1), beta(2), beta(3), beta(4)), guessbeta, lowerbound, upperbound, options);
        
        % compute the x and y values for plotting
        % fit the cumulative Weibull
        fit2.y = cumWB(fit1.x, fit2.mean.beta(1), fit2.mean.beta(2), fit2.mean.beta(3), fit2.mean.beta(4));
        semilogx(fit1.x*100, fit2.y, '-k', 'LineWidth', 2);
        
        % add param info
        text(20, .5, gca, sprintf('threshold %.2f %% \n slope %.2f \n guessrate %.2f %% \n lapse rate %.2f %%', ...
            fit2.mean.beta(1)*100, fit2.mean.beta(2), fit2.mean.beta(3)*100, fit2.mean.beta(4)*100)); hold on;
        
        % ----------- PARAMETRIC BOOTSTRAP FOR WEIBULL PARAMETER CIs -------
        
        switch bootstrapping
            case 'bootstrap',
                for iboot = 1:nboot,
                    
                    % get the percent correct based on the weibull fit
                    prob = cumWB(datapoints1.all.difficulty, fit2.mean.beta(1), fit2.mean.beta(2), fit2.mean.beta(3), fit2.mean.beta(4));
                    % create a mock response set
                    % this is the source of randomness in the bootstrap
                    boot_resp = binornd(ones(1,size(datapoints2.all.response,2)), prob(1:size(datapoints2.all.response,2)));
                    
                    % change this into parameters for weibull fit
                    for idiff = 1:length(datapoints1.mean.difficulty),
                        thisdiff = datapoints1.mean.difficulty(idiff);
                        % parametric bootstrap for weibull parameters
                        boot_corr(idiff) = nanmean(length(find(datapoints1.all.difficulty == thisdiff & boot_resp == 1))/length(find(datapoints1.all.difficulty==thisdiff)));
                        
                    end % coh
                    
                    % fit the weibull function to this set of betas
                    [fit2.boot.beta(:,iboot), fit2.boot.fval(iboot)] = fminsearchbnd(@(beta) LL_cumWB(datapoints1.mean.difficulty, ...
                        boot_corr, beta(1), beta(2), beta(3), beta(4)), guessbeta, lowerbound, upperbound, options);
                    
                end % boot
                
            case 'nobootstrap'
                %  do nothing
        end
end

% plot this
semilogx(datapoints1.mean.difficulty*100, datapoints2.mean.correctedy, 'ro', 'MarkerSize',8,  'MarkerFaceColor', 'r');

switch bootstrapping,
    case 'bootstrap',
        % and the confidence interval errorbars
        datapoints2.bootstrap.CI.lower = prctile(datapoints2.bootstrap.correctedy', fit1.CI(1));
        datapoints2.bootstrap.CI.upper = prctile(datapoints2.bootstrap.correctedy', fit1.CI(2));
        errorbar(datapoints1.mean.difficulty*100, datapoints2.mean.correctedy, ...
            datapoints2.mean.correctedy - datapoints2.bootstrap.CI.lower, datapoints2.bootstrap.CI.upper - datapoints2.mean.correctedy, 'r.');
        
        % just take the percentiles
        fit2.boot.CI.lower = prctile(fit2.boot.beta(1,:)', fit1.CI(1));
        fit2.boot.CI.upper = prctile(fit2.boot.beta(1,:)', fit1.CI(2));
end

%plot
title(sprintf('D''-corrected performance'));
xlabel('Motion coherence from baseline (%)'); ylabel('Corrected accuracy (% correct)');
ylim([.39 1.01]);
set(gca, 'XTick', datapoints1.mean.difficulty*100);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSYCHOMETRIC FUNCTION NR3 - RANGE OF COHERENCES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(223);

% reshape
datapoints3.mean.coherence           = unique(datapoints2.all.coherence);
datapoints3.all.response             = datapoints2.all.response;
datapoints3.all.response(datapoints3.all.response == 0) = -1; %change back so I can take the mean
for icoh = 1:length(datapoints3.mean.coherence),
    datapoints3.mean.response(icoh)             = nanmean(datapoints3.all.response(datapoints2.all.coherence==datapoints3.mean.coherence(icoh)));
end

% ------- BOOTSTRAP DATAPOINTS ---------------------------------------------------------------------------
switch bootstrapping,
    case 'bootstrap',
        disp('bootstrapping...');
        
        % preallocate
        datapoints3.bootstrap.response = nan(length(datapoints3.mean.coherence), nboot);
        
        for iboot = 1:nboot,
            
            % for each coherence level, calculate the proportion of 'stronger'
            % responses
            for icoh = 1:length(datapoints3.mean.coherence),
                
                % get a random sample from all the datapoints with replacement -
                % only for this level of coherence
                [bootstrappeddata] = datasample(datapoints3.all.response(datapoints2.all.coherence == datapoints3.mean.coherence(icoh)), ...
                    length(datapoints2.all.coherence(datapoints2.all.coherence == datapoints3.mean.coherence(icoh))), 'Replace', true);
                
                datapoints3.bootstrap.response(icoh, iboot) = mean(bootstrappeddata);
            end
        end
        % normalize to range 0-1
        datapoints3.bootstrap.response  = (datapoints3.bootstrap.response+1)/2;
        
end

% normalize to range 0-1
datapoints3.mean.response       = (datapoints3.mean.response+1)/2;

% ------- FIT LOGISTIC FUNCTION  ---------------------------------------------------------------------------
switch fitting,
    case 'fit',
        
        logistic = @(coh, threshold, slope, lower, upper) lower+(1-lower-upper)*(1./(1+exp(-1*(slope).*(coh-threshold))));
        
        LL_logistic = @(intensities, results, threshold, slope, lower, upper) -sum(results.*log(logistic(intensities, ...
            threshold, slope, lower, upper)) + (1-results).*log(1-logistic(intensities, threshold, slope, lower, upper)));
        
        % initial values for fminsearch
        guessbeta(1) = .70; %PSE
        guessbeta(2) = 10; %sensitivity, slope
        guessbeta(3) = 0; %lower rate
        guessbeta(4) = 0; % upper rate (lapse rate)
        
        % lower bound
        lowerbound(1)   = 0;
        lowerbound(2)   = 0;
        lowerbound(3)   = 0;
        lowerbound(4)   = 0;
        
        % upper bound
        upperbound(1)   = +inf;
        upperbound(2)   = +inf;
        upperbound(3)   = .6;
        upperbound(4)   = .1;
        
        % find optimal mean values for beta using fminsearch
        [fit3.mean.beta, fit3.mean.fval] = fminsearchbnd(@(beta) LL_logistic(datapoints3.mean.coherence, ...
            datapoints3.mean.response, beta(1), beta(2), beta(3), beta(4)), guessbeta, lowerbound, upperbound, options);
        
        % ----------- PARAMETRIC BOOTSTRAP FOR WEIBULL PARAMETER CIs -------
        
        switch bootstrapping
            case 'bootstrap',
                for iboot = 1:nboot,
                    
                    % get the percent correct based on the weibull fit
                    prob = cumWB(datapoints2.all.coherence, fit3.mean.beta(1), fit3.mean.beta(2), fit3.mean.beta(3), fit3.mean.beta(4));
                    % create a mock response set
                    % this is the source of randomness in the bootstrap
                    boot_resp = binornd(ones(1,size(datapoints3.all.response,2)), prob(1:size(datapoints3.all.response,2)));
                    
                    % change this into parameters for weibull fit
                    for icoh = 1:length(datapoints3.mean.coherence),
                        thiscoh = datapoints3.mean.coherence(icoh);
                        % parametric bootstrap for weibull parameters
                        boot_corr(icoh) = nanmean(length(find(datapoints2.all.coherence == thiscoh & boot_resp == 1))/length(find(datapoints2.all.coherence == thiscoh)));
                    end % coh
                    
                    % fit the logistic function to this set of betas
                    [fit3.boot.beta(:,iboot), fit3.boot.fval(iboot)] = fminsearchbnd(@(beta) LL_logistic(datapoints3.mean.coherence, ...
                        boot_corr, beta(1), beta(2), beta(3), beta(4)), guessbeta, lowerbound, upperbound, options);
                    
                end % boot
                
            case 'nobootstrap'
                %  do nothing
        end
        
        % ------- PLOT  ---------------------------------------------------------------------------
        
        % compute the x and y values for plotting
        fit3.x = linspace(min(datapoints3.mean.coherence),max(datapoints3.mean.coherence));
        % fit the cumulative Weibull
        fit3.y = logistic(fit3.x, fit3.mean.beta(1), fit3.mean.beta(2), fit3.mean.beta(3), fit3.mean.beta(4));
        
        %plot
        plot(fit3.x*100, fit3.y, '-k', 'LineWidth', 2); hold on;
end

% plot the actual datapoints on top
plot(datapoints3.mean.coherence*100, datapoints3.mean.response, 'go', 'MarkerSize',10, 'MarkerFaceColor', 'g'); hold on;

switch bootstrapping,
    case 'bootstrap'
        % and the confidence interval errorbars
        datapoints3.bootstrap.CI.lower = prctile(datapoints3.bootstrap.response', fit1.CI(1));
        datapoints3.bootstrap.CI.upper = prctile(datapoints3.bootstrap.response', fit1.CI(2));
        errorbar(datapoints3.mean.coherence*100, datapoints3.mean.response, ...
            datapoints3.mean.response - datapoints3.bootstrap.CI.lower, datapoints3.bootstrap.CI.upper - datapoints3.mean.response, 'g.');
        
        % just take the percentiles
        fit3.boot.CI.lower = prctile(fit3.boot.beta(2,:)', fit1.CI(1));
        fit3.boot.CI.upper = prctile(fit3.boot.beta(2,:)', fit1.CI(2));        
end

title(sprintf('Full response'));
xlabel('Coherence (%)'); ylabel('% responses stronger motion');
ylim([-.01 1.01]); %xlim([38 102]);
%set(gca, 'XTick', datapoints3.mean.coherence*100, 'XTickLabel', round(datapoints3.mean.coherence*100));

switch fitting,
    case 'fit'
        % add param info
        text(80, .2, sprintf('PSE %.2f %% \n sensitivity %.2f \n guessrate %.2f %% \n lapse rate %.2f %%', ...
            fit3.mean.beta(1)*100, fit3.mean.beta(2), fit3.mean.beta(3)*100, fit3.mean.beta(4)*100));
end
% ------- GENERAL NICE THINGS FOR THE WHOLE PLOT  ---------------------------------------------------------------------------

suptitle(sprintf('Subject %d, session %d', setup.participant, setup.session));

% % make it all look a bit better
set(findall(gca,'type','text'),'fontWeight', 'bold', 'FontSize', 14, 'FontName', 'Helvetica');

% make sure the output doesn't go awry
switch fitting,
    case 'nofit',
        fit1 = [];
        fit2 = [];
        fit3 = [];
end

% save the fig
%set(gcf,'PaperOrientation','landscape');
%set(gcf,'PaperUnits','normalized');
%set(gcf,'PaperPosition', [0 0 1 1]);
%saveas(gcf, sprintf('PsychFuncFit_P%d_session%d.eps', setup.participant, setup.session), 'epsc');

%save the data
%save(sprintf('PsychFuncFit_P%d_session%d.mat', setup.participant, setup.session), 'fit1', 'fit2', 'fit3');


end