function mediationAnalysis()
% does a mediation analysis and the Sobell-Aroian test to see if pupil
% diameter and/or RT mediate the effect of uncertainty on switching
% behaviour

global mypath;

% do this in lavaan (R) or manually in Matlab
% when using probit regression, these give the same results
whichLanguage = 'R';
% whichLanguage = 'Matlab';

% ============================================ %
% write files
% ============================================ %

for sj = 1:27,
  %  if ~exist(sprintf('%s/Data/CSV/SEMdata_sj%02d.csv', mypath, sj), 'file'),
        
        disp(sj);
        clearvars -except mypath grandavg sj whichLanguage;
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
        data.rt                 = zscore(log(data.rt + 0.1));
        
        % =================================================== %
        % step 2b: write 2 csv for loading into lavaan
        % =================================================== %
        
        writetable(data, sprintf('%s/Data/CSV/SEMdata_sj%02d.csv', mypath, sj));
   % end
end

% =================================================== %
% test correlation between uncertainty and pupil
% =================================================== %

for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/SEMdata_sj%02d.csv', mypath, sj));
    [grandavg.rho(sj), grandavg.pval(sj)] = corr(data.uncertainty, data.decision_pupil);
end

fprintf('r = %.3f, range %.3f to %.3f, significant in %d out of %d participants \n', ...
    mean(grandavg.rho), min(grandavg.rho), max(grandavg.rho), length(find(grandavg.pval < 0.05)), 27);

% =================================================== %
% run actual mediation
% =================================================== %

switch whichLanguage,
    case 'Matlab'
        
        for sj = 1:27,
            data = readtable(sprintf('%s/Data/CSV/SEMdata_sj%02d.csv', mypath, sj));
            
            % if we want to compare to R, better to use the same link function
            % (lavaan doesn't support logistic regression)
            whichFit = 'probit';
            
            % =================================================== %
            % step 3: run all the regression
            % =================================================== %
            
            % eq 1. direct effect of X on Y, binary outcome
            mdl = fitglm(data, 'repeat ~ 1 + stimrepeat + uncertainty', ...
                'distr', 'binomial', 'link', whichFit);
            grandavg.coef1(sj, :)   = mdl.Coefficients.Estimate;
            grandavg.std1(sj, :)    = mdl.Coefficients.SE;
            grandavg.eq1Names       = mdl.CoefficientNames; % these wont be in the order of input
            
            % eq 2. direct effect of M on Y, linear
            mdl = fitglm(data, 'decision_pupil ~ 1 + uncertainty');
            grandavg.coef2(sj, :)   = mdl.Coefficients.Estimate;
            grandavg.std2(sj, :)    = mdl.Coefficients.SE;
            grandavg.eq2Names       = mdl.CoefficientNames;
            
            % eq 3. effect of X and M together on Y, binary outcome
            mdl = fitglm(data, 'repeat ~ 1 + stimrepeat + uncertainty + decision_pupil', ...
                'distr', 'binomial', 'link', whichFit);
            grandavg.coef3(sj, :)   = mdl.Coefficients.Estimate;
            grandavg.std3(sj, :)    = mdl.Coefficients.SE;
            grandavg.eq3Names       = mdl.CoefficientNames;
            
            % extract the relevant information
            grandavg.c(sj)   = grandavg.coef1(sj, strcmp('uncertainty',     grandavg.eq1Names))'; % without b
            grandavg.a(sj)   = grandavg.coef2(sj, strcmp('uncertainty',     grandavg.eq2Names))';
            grandavg.c1(sj)  = grandavg.coef3(sj, strcmp('uncertainty',     grandavg.eq3Names))'; % with b
            grandavg.b(sj)   = grandavg.coef3(sj, strcmp('decision_pupil',  grandavg.eq3Names))';
            
            % =================================================== %
            % compute standard deviation of all parameters
            % see % http://davidakenny.net/doc/dichmed.pdf for these derivations
            % =================================================== %
            
            grandavg.stdX(sj)       = std(data.uncertainty);
            grandavg.stdM(sj)       = std(data.decision_pupil);
            
            % if we want to compare to R, better to use the same link function
            % (lavaan doesn't support logistic regression)
            % see % http://davidakenny.net/doc/dichmed.pdf
            switch whichFit
                case 'probit'
                    factorScale = 1;
                case 'logit'
                    factorScale = pi^2/3;
            end
            
            % Var(Y') = c^2 * Var(X) + p^2/3
            grandavg.stdY1(sj)      = sqrt(grandavg.c(sj)^2 * grandavg.stdX(sj)^2 + factorScale)';
            
            % Var(M') = a^2 * Var(X) + p^2/3
            % we won't use this one actually, since a is estimated by OLS
            grandavg.stdM1(sj)      = sqrt(grandavg.a(sj)^2 * grandavg.stdX(sj)^2 + factorScale)' ;
            
            % Var(Y") = c'^2 * Var(X) + b^2 * Var(M) + 2*b*c'*Cov(X,M) + p^2/3
            % the covariance is scalar, so take the Pearson correlation
            grandavg.stdY2(sj)      = sqrt(grandavg.c1(sj)^2 *grandavg.stdX(sj)^2 ...
                + grandavg.b(sj)^2 * grandavg.stdM(sj)^2 ...
                + 2 * grandavg.b(sj) * grandavg.c1(sj) * corr(data.uncertainty(:), data.decision_pupil(:)) ...
                +  factorScale');
        end
        
        % =================================================== %
        % compute the corrected coefficients by multiplying by the standard
        % deviation of the predictor variable, then dividing by the standard
        % deviation of the outcome variable.
        % =================================================== %
        
        dat     = table;
        dat.a   = grandavg.a'; % a is an OLS coefficient and does not need rescaling
        dat.b   = transpose(grandavg.b .* grandavg.stdM ./ grandavg.stdY2);
        dat.c   = transpose(grandavg.c .* grandavg.stdX ./ grandavg.stdY1);
        dat.c1  = transpose(grandavg.c1 .* grandavg.stdX ./ grandavg.stdY2);
        
        % also save the paths, direct and indirect
        dat.indirect  = dat.a .* dat.b;
        dat.indirect2 = dat.c - dat.c1;
        dat.direct    = dat.c;
        dat.total     = dat.c1 + dat.a .* dat.b;
        
        % check that the indirect paths computed in those two ways are
        % identical up to numerical tolerance
        assert(all(abs(dat.indirect - dat.indirect2) < 0.001), 'indirect paths not correctly computed');
        writetable(dat, sprintf('%s/Data/CSV/SEMdata_Matlab.csv', mypath));
        
    case 'R'
        % call R from the command line, point to absolute path
        % an output file will be saved from R
        system('/Library/Frameworks/R.framework/Resources/bin/R < MediationLavaan.R --no-save');
        dat = readtable(sprintf('%s/Data/CSV/SEMdata_lavaan.csv', mypath));
        
        %plotBetasSwarm(dat{:, 2:end})
        %set(gca, 'xtick', 1:11, 'xticklabel', dat.Properties.VariableNames(2:end))
end

% =================================================== %
% visualise
% =================================================== %

clf; subplot(4,4,1);
plotBetasSwarm([dat.a dat.b dat.c1]);
set(gca, 'xtick', 1:5, 'xticklabel', {'a', 'b', 'c'''});

subplot(4,8,4);
plotBetasSwarm([dat.indirect]);
set(gca, 'xtick', 1:2, 'xticklabel', {'a1*b1', 'a2*b2'});

print(gcf, '-dpdf', sprintf('%s/Figures/mediationAnalysis.pdf', mypath));

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

[h(1), pval(1)] = ttest(dat.c);

% Step 2: Show that the causal variable is correlated with the mediator.
% Use M as the criterion variable in the regression equation and X as a predictor
% (estimate and test path a).  This step essentially involves treating the mediator
% as if it were an outcome variable.

[h(2), pval(2)] = ttest(dat.a);

% Step 3:  Show that the mediator affects the outcome variable.  Use Y as the
% criterion variable in a regression equation and X and M as predictors
% (estimate and test path b).  It is not sufficient just to correlate the
% mediator with the outcome because the mediator and the outcome may be
% correlated because they are both caused by the causal variable X.  Thus,
% the causal variable must be controlled in establishing the effect of the mediator on the outcome.

[h(3), pval(3)] = ttest(dat.b);

% Step 4:  To establish that M completely mediates the X-Y relationship,
% the effect of X on Y controlling for M (path c') should be zero (see discussion
% below on significance testing). The effects in both Steps 3 and 4 are
% estimated in the same equation.

% test if c1 is closer to zero than c
[h(4), pval(4)] = ttest(dat.c, dat.c1);

% test if all these conditions are met!
if all(pval < 0.05),
    disp('all causal steps checked off!');
else
    disp('mediation not supported');
end

assert(all(pval < 0.05), 'mediation not supported');
disp('done!');

% =================================================== %
% use model comparison to see if the model with mediation
% does better than the model without
% =================================================== %


end
