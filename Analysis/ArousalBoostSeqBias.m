%% JWs data with confidence, arousal boost and RT - predict sequential effects

clear; clc; close;
global mypath;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';
warning('error', 'stats:glmfit:IterationLimit');

alldata         = readtable(sprintf('%s/Data/CSV/arousalboost.csv', mypath));
alldata(1, :)   = []; % remove idx;
subjects        = unique(alldata.subj_idx);

for sj = subjects',
    disp(sj);
    % get data for this subject
    data = alldata((alldata.subj_idx == sj), :);
    
    % make previous trial regressors
    data.prevRT         = circshift(zscore(log(data.rt + 0.1)), 1);
    
    % effect code binary variables
    data.prevResp       = circshift(data.response, 1);
    data.prevResp(data.prevResp == 0) = -1;
    data.prevStim       = circshift(data.stimulus, 1);
    data.prevStim(data.prevStim == 0) = -1;
    data.prevBoost      = circshift(data.split, 1);
    data.prevBoost(data.prevBoost == 0) = -1;
    data.prevConf       = circshift(data.confidence, 1);
    data.prevConf(data.prevConf == 0) = -1;
    
    % design matrix
    % eq 1. direct effect of X on Y, logistic
    try
        mdl = fitglm(data, ['response ~ 1 + stimulus + prevResp + prevStim' ...
            '+ prevResp*prevRT + prevStim*prevRT + prevRT'...
            '+ prevResp*prevBoost + prevStim*prevBoost + prevBoost' ...
            '+ prevResp*prevConf + prevStim*prevConf + prevConf'], ...
            'distr', 'binomial', 'link', 'logit');
        
        grandavg.coef(sj+1, :) = mdl.Coefficients.Estimate;
        grandavg.eqNames       = mdl.CoefficientNames; % these wont be in the order of input
    catch
        % remove confidence from the model
        mdl = fitglm(data, ['response ~ 1 + stimulus + prevResp + prevStim' ...
            '+ prevResp*prevRT + prevStim*prevRT + prevRT'...
            '+ prevResp*prevBoost + prevStim*prevBoost + prevBoost'], ...
            'distr', 'binomial', 'link', 'logit');
        grandavg.coef(sj+1, [1:6 8:11]) = mdl.Coefficients.Estimate;
        grandavg.coef(sj+1, [7 12 13]) = NaN;
    end
 
end

% table
grandavg.eqNames = regexprep(grandavg.eqNames, '(', '');
grandavg.eqNames = regexprep(grandavg.eqNames, ')', '');
grandavg.eqNames = regexprep(grandavg.eqNames, ':', 'X');
dat = array2table(grandavg.coef, 'variablenames', grandavg.eqNames);

%%
clf;
% betas
subplot(4,4,[1:3]);
plotBetasSwarm(grandavg.coef);
set(gca, 'xtick', 1:length(grandavg.eqNames), ...
    'xticklabel', grandavg.eqNames, 'xticklabelrotation', -30);

% decision strategies
subplot(4,4,4);
hold on;
plot([-1 1], [-1 1], 'color', 'k', 'linewidth', 0.5);
plot([-1 1], [1 -1], 'color', 'k', 'linewidth', 0.5);
scatter(dat.prevResp, dat.prevStim, 10, 'filled');
xlabel('Choice weight'); ylabel('Stimulus weight');
box on; axis square; axis tight;

print(gcf, '-dpdf', sprintf('%s/Figures/arousalBoost.pdf', mypath));
