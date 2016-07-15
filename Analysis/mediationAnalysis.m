function mediationAnalysis()
% does a mediation analysis and the Sobell-Aroian test to see if pupil
% diameter and/or RT mediate the effect of uncertainty on switching
% behaviour

% ============================================ %
% check if uncertainty behaves as we expect
% ============================================ %

%nbins = 6;
%subplot(221); psychFuncShift_Bias('uncertainty', nbins, 1); title('Correct');
%subplot(222); psychFuncShift_Bias('uncertainty', nbins, 1); title('Error');

global mypath;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

if 0,
    for sj = 1:27,
        disp(sj);
        clearvars -except mypath grandavg sj;
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        data = data(find(data.sessionnr > 1), :);
        
        % ============================================ %
        % step 1: compute model-based uncertainty
        % ============================================ %
        
        % fit a probit slope so we can get the sigma
        b       = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');
        sigma   = 1/b(2);  % standard deviation at these values is the inverse!
        bound   = -b(1); % the bound is negative if people say
        
        % plot
        % yfit  = glmval(b,data.motionstrength,'probit');
        % plot(data.motionstrength, yfit, '.');
        % grid on; xlim([-6 6]);ylim([0 1]);
        
        % for each trial, compute the average level of uncertainty
        data.uncertainty = arrayfun(@simulateUncertainty, abs(data.motionstrength), ...
            data.correct, sigma*ones(length(data.correct), 1), bound*ones(length(data.correct), 1));
        
        % =================================================== %
        % step 2a: prepare data for repetition model
        % =================================================== %
        
        repetition          = (diff(data.resp) == 0);
        data.repeat         = [repetition; NaN];
        
        % measure of how much the stimulus tells the subject to repeat from
        % their last choice to this one
        % important: the logistic slope for this measure of stimulus repetition
        % is steeper than for the one below!
        stimrepeat          = data.stim - circshift(data.resp, -1);
        stimrepeat          = double(stimrepeat == 0);
        stimrepeat(find(stimrepeat == 0)) = -1;
        stimrepeat          = stimrepeat .* circshift(abs(data.motionstrength), -1);
        data.stimrepeat     = stimrepeat;
        
        % don't use trials that do not have the subsequent trial following them
        % (and those at the end of each block)
        trlDif                  = [diff(data.trialnr); 0];
        removeTrls              = false(size(trlDif));
        removeTrls(trlDif < 1)  = true;
        removeTrls(trlDif > 1)  = true;
        data.resp(removeTrls) 	= NaN; % avoid modelling these
        
        % pupil scores have already been normalized
        % zscore the uncertainty values, these are NOT normally distributed!
        data.decision_pupil     = zscore(data.decision_pupil);
        data.uncertainty        = zscore(data.uncertainty);
        
        % =================================================== %
        % step 2b: write 2 csv for loading into lavaan
        % =================================================== %
        
        writetable(data, sprintf('%s/Data/CSV/SEMdata_sj%02d.csv', mypath, sj));
    end
end

% choose
whichFit = 'probit';

for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/SEMdata_sj%02d.csv', mypath, sj));
    
    % =================================================== %
    % step 3: run all the regression
    % =================================================== %
    
    try
    % eq 1. direct effect of X on Y, logistic
    mdl = fitglm(data, 'repeat ~ 1 + stimrepeat + uncertainty', ...
        'distr', 'binomial', 'link', whichFit);
    catch
        assert(1==0);
    end
    grandavg.coef1(sj, :)   = mdl.Coefficients.Estimate;
    grandavg.std1(sj, :)    = mdl.Coefficients.SE;
    grandavg.eq1Names       = mdl.CoefficientNames; % these wont be in the order of input
    
    % eq 2. direct effect of M on Y, linear
    mdl = fitglm(data, 'decision_pupil ~ 1 + uncertainty');
    grandavg.coef2(sj, :)   = mdl.Coefficients.Estimate;
    grandavg.std2(sj, :)    = mdl.Coefficients.SE;
    grandavg.eq2Names       = mdl.CoefficientNames;
    
    % eq 3. effect of X and M together on Y, logistic
    mdl = fitglm(data, 'repeat ~ 1 + stimrepeat + uncertainty + decision_pupil', ...
        'distr', 'binomial', 'link', whichFit);
    grandavg.coef3(sj, :)   = mdl.Coefficients.Estimate;
    grandavg.std3(sj, :)    = mdl.Coefficients.SE;
    grandavg.eq3Names       = mdl.CoefficientNames;
    
    % extract the relevant information
    grandavg.c(sj)   = grandavg.coef1(sj, strcmp('uncertainty',     grandavg.eq1Names)); % without b
    grandavg.a(sj)   = grandavg.coef2(sj, strcmp('uncertainty',     grandavg.eq2Names));
    grandavg.c1(sj)  = grandavg.coef3(sj, strcmp('uncertainty',     grandavg.eq3Names)); % with b
    grandavg.b(sj)   = grandavg.coef3(sj, strcmp('decision_pupil',  grandavg.eq3Names));
    
    % to compute the Sobel test, from
    
    grandavg.stdX(sj)       = std(data.uncertainty);
    grandavg.stdM(sj)       = std(data.decision_pupil);
    
    % see % http://davidakenny.net/doc/dichmed.pdf for these derivations
    % note: use 1 instead of pi^2/3 when model is probit
    
    switch whichFit
        case 'probit'
            factorScale = 1;
        case 'logit'
            factorScale = pi^2/3;
    end
    
    % Var(Y') = c^2 * Var(X) + p^2/3
    grandavg.stdY1(sj)      = sqrt(grandavg.c(sj)^2 * grandavg.stdX(sj)^2 + factorScale);
    
    % Var(M') = a^2 * Var(X) + p^2/3
    % we won't use this one actually, since a is estimated by OLS
    grandavg.stdM1(sj)      = sqrt(grandavg.a(sj)^2 * grandavg.stdX(sj)^2 + factorScale) ;
    
    % Var(Y") = c'^2 * Var(X) + b^2 * Var(M) + 2*b*c'*Cov(X,M) + p^2/3
    % the covariance is scalar, so take the Pearson correlation
    grandavg.stdY2(sj)      = sqrt(grandavg.c1(sj)^2 *grandavg.stdX(sj)^2 ...
        + grandavg.b(sj)^2 * grandavg.stdM(sj)^2 ...
        + 2 * grandavg.b(sj) * grandavg.c1(sj) * corr(data.uncertainty(:), data.decision_pupil(:)) ...
        +  factorScale);
end

% visualise
figure; subplot(4,4,1);
plotBetasSwarm([grandavg.a' grandavg.b']);
set(gca, 'xtick', 1:2, 'xticklabel', {'a', 'b'});
subplot(4,4,2);
plotBetasSwarm([grandavg.c' grandavg.c1']);
set(gca, 'xtick', 1:2, 'xticklabel', {'c', 'c'''});

% =================================================== %
% step 4: causal steps approach, Baron and Kenny 1986
% from http://davidakenny.net/cm/mediate.htm
% =================================================== %
disp('testing the causal steps approach...');

% also get the standard error on each

% Step 1:  Show that the causal variable is correlated with the outcome.
% Use Y as the criterion variable in a regression equation and X as a predictor
% (estimate and test path c in the above figure). This step establishes that
% there is an effect that may be mediated.

[h(1), pval(1)] = ttest(grandavg.c);

% Step 2: Show that the causal variable is correlated with the mediator.
% Use M as the criterion variable in the regression equation and X as a predictor
% (estimate and test path a).  This step essentially involves treating the mediator
% as if it were an outcome variable.

[h(2), pval(2)] = ttest(grandavg.a);

% Step 3:  Show that the mediator affects the outcome variable.  Use Y as the
% criterion variable in a regression equation and X and M as predictors
% (estimate and test path b).  It is not sufficient just to correlate the
% mediator with the outcome because the mediator and the outcome may be
% correlated because they are both caused by the causal variable X.  Thus,
% the causal variable must be controlled in establishing the effect of the mediator on the outcome.

[h(3), pval(3)] = ttest(grandavg.b);

% Step 4:  To establish that M completely mediates the X-Y relationship,
% the effect of X on Y controlling for M (path c') should be zero (see discussion
% below on significance testing). The effects in both Steps 3 and 4 are
% estimated in the same equation.

% test if c1 is closer to zero than c
[h(4), pval(4)] = ttest(grandavg.c, grandavg.c1);

% test if all these conditions are met!
if all(pval < 0.05),
    disp('all causal steps checked off!');
else
    disp('mediation not supported');
end

assert(all(pval < 0.05), 'mediation not supported');
disp('done!');

% =================================================== %
% step 5: measure proportion explained by the indirect effect
% from http://davidakenny.net/cm/mediate.htm
% =================================================== %

% compute the corrected coefficients by multiplying by the standard
% deviation of the predictor variable, then dividing by the standard
% deviation of the outcome variable.

grandavg.a_corr  = grandavg.a;
% a is an OLS coefficient and does not need rescaling
grandavg.b_corr  = grandavg.b .* grandavg.stdM ./ grandavg.stdY2;
grandavg.c_corr  = grandavg.c .* grandavg.stdX ./ grandavg.stdY1;
grandavg.c1_corr = grandavg.c1 .* grandavg.stdX ./ grandavg.stdY2;

fprintf('Test: mean ab+c1 = %.3f, c = %.3f, difference %.3f \n', ...
    mean(grandavg.a_corr .* grandavg.b_corr + grandavg.c1_corr), ...
    mean(grandavg.c_corr), mean(grandavg.c_corr - (grandavg.a_corr .* grandavg.b_corr + grandavg.c1_corr)));

% also the standard errors of these rescaled coefficients for inference
grandavg.a_corr_se = grandavg.std2(:, strcmp('uncertainty', grandavg.eq2Names));
grandavg.b_corr_se = grandavg.std2(:, strcmp('uncertainty', grandavg.eq2Names));

% % corrected plot
% subplot(4,4,3);
% plotBetasSwarm([grandavg.a' grandavg.b_corr']);
% set(gca, 'xtick', 1:2, 'xticklabel', {'a', 'b'});
% title('Normalized');
% subplot(4,4,4);
% plotBetasSwarm([grandavg.c_corr' grandavg.c1_corr']);
% set(gca, 'xtick', 1:2, 'xticklabel', {'c', 'c'''});
% title('Normalized');

subplot(4,4,3);
histogram([grandavg.c - grandavg.c1], 15);
xlabel('c - c'''); ylabel('count'); box off;

% show indirect effect
indirect = grandavg.a_corr .* grandavg.b_corr;
indirect2 = grandavg.c_corr - grandavg.c1_corr;

subplot(4,4,4);
plotBetasSwarm(indirect2');
xlabel('indirect effect');
set(gca, 'xtick', []);

print(gcf, '-dpdf', sprintf('%s/Figures/mediationAnalysis.pdf', mypath));

% save
savefast(sprintf('%s/Data/GrandAverage/mediationModel.mat', mypath), 'grandavg');

% =================================================== %
% use the Sobell test to see if the indirect path is significant
% from http://quantpsy.org/sobel/sobel.htm
% =================================================== %
if 0,
sobel.z = grandavg.a * grandavg.b / ...
    sqrt(grandavg.b .^2 * sa.^2 + a .^2 * sb .^2); % original Sobel test
sobel.z = grandavg.a * grandavg.b / sqrt(b.^2*sa.^2 + a.^2*sb.^2 + sa.^2*sb.^2); % Aroian test
sobel.p = normcdf(-abs(sobel.z), 0, 1); % get the p-value for this z-value;
end

% =================================================== %
% compare Lavaan with my own mediation parameters
% =================================================== %

clearvars -except mypath
load(sprintf('%s/Data/GrandAverage/mediationModel.mat', mypath));
rdat = readtable(sprintf('%s/Data/CSV/SEMdata_lavaan.csv', mypath));

figure;
% plot a scatterHistDiff for each parameter
vars = {'a', 'b', 'c', 'c1'};
for v = 1:length(vars),
    subplot(4,4,v);
    %      scatterHistDiff(rdat.(vars{v})', grandavg.(vars{v}));
    scatterHistDiff(rdat.(vars{v})', grandavg.([vars{v} '_corr']));
    title(vars{v});
end

end
