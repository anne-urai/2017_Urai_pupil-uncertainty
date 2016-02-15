%% goodness of fit for logistic models

clear;
subjects = 1:27; clc
for sj = unique(subjects),
    
    data = readtable(sprintf('~/Data/pupilUncertainty/Data/CSV/2ifc_data_sj%02d.csv', sj));
    
    data.resp(data.resp==-1) = 0;
    
    % first, normalize both;
    data.cohstim            = zscore(data.coherence .* data.stim);
    b                       = regress(data.cohstim, [ones(size(data.motionstrength)) data.motionstrength]);
    data.motionstrength     = data.motionstrength * b(2) + b(1);
    
    % first, fit psychfunc on coherence
    mdlC                    = fitglm(data, 'resp ~ cohstim', 'Distribution','binomial');
    grandavg.bC(sj, :)      = mdlC.Coefficients.Estimate;
    grandavg.logLC(sj)      = mdlC.LogLikelihood;
    grandavg.bicC(sj)       = -2*(grandavg.logLC(sj))+2 .* log(height(data));
    
    % now motionenergy
    mdlM                    = fitglm(data, 'resp ~ motionstrength', 'Distribution','binomial');
    grandavg.bM(sj, :)      = mdlM.Coefficients.Estimate;
    grandavg.logLM(sj)      = mdlM.LogLikelihood;
    grandavg.bicM(sj)       = -2*(grandavg.logLM(sj))+2 .* log(height(data));
    
    subplot(5,6,sj);
    x = -3:0.01:3;
    plot(x, mdlC.predict(x))
    
end

% then compare those two
[h, p] = ttest(grandavg.bicM, grandavg.bicC);
BICdiff = grandavg.bicM - grandavg.bicC;
% on average, the BIC for M is 7 higher than the BIC for C

[h, p] = ttest(grandavg.logLM, grandavg.logLC);
logLdiff = grandavg.logLM - grandavg.logLC;
% on average, the logL for M is 3.5 smaller than logC for C

beta = [grandavg.bM(:, 2) grandavg.bC(:, 2)];
[h, p] = ttest(beta(:, 1), beta(:, 2));
% the beta for M is -0.027 smaller than the beta for C

%% another way: within one coherence level, split by filtered motionenergy
