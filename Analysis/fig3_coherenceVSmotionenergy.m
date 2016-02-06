%% goodness of fit for logistic models

clear;
subjects = 1:27; clc
for sj = unique(subjects),
    
    data = readtable(sprintf('~/Data/pupilUncertainty/Data/CSV/2ifc_data_sj%02d.csv', sj));
    
    data.resp(data.resp==-1) = 0;
    
    % first, fit psychfunc on coherence
    data.cohstim            = zscore(data.coherence .* data.stim);
    mdlC                    = fitglm(data, 'resp ~ cohstim', 'Distribution','binomial');
    grandavg.bC(sj, :)      = mdlC.Coefficients.Estimate;
    grandavg.logLC(sj)      = mdlC.LogLikelihood;
    grandavg.bicC(sj)       = -2*(grandavg.logLC(sj))+2 .* log(height(data));
    
    % now motionenergy
    data.motionstrength     = zscore(data.motionstrength);
    mdlM                    = fitglm(data, 'resp ~ motionstrength', 'Distribution','binomial');
    grandavg.bM(sj, :)      = mdlM.Coefficients.Estimate;
    grandavg.logLM(sj)      = mdlM.LogLikelihood;
    grandavg.bicM(sj)       = -2*(grandavg.logLM(sj))+2 .* log(height(data));
    
end

%% another way: within one coherence level, split by filtered motionenergy
