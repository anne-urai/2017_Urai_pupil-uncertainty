%% redo the whole thing in one big mixed effects model
% Moscatelli et al. JoV

clear all; clc; close all;

data = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data2_allsj.csv'));

% create a unique id for each sj/session combination
newsession = nan(size(data.subjnr));
cnt = 1;
for sj = unique(data.subjnr)',
    sessions = unique(data.sessionnr(sj==data.subjnr));
    for session = sessions',
        newsession(data.subjnr==sj & data.sessionnr == session) = (cnt);
        cnt = cnt + 1;
    end
end

%newsession2 = num2str(newsession);
data.newsession = newsession;

% normalize the pupil response for each subject
% this way, we can interpret the coefficient of the pupil fixed effect in
% units of standard deviation
for sj = unique(data.newsession)',
    data.decision_pupil(sj==data.newsession) = zscore(data.decision_pupil(sj==data.newsession));
    data.motionstrength(sj==data.newsession) = zscore(data.motionstrength(sj==data.newsession));
    data.rt(sj==data.newsession) = zscore(log(data.rt(sj==data.newsession)));
end

% history variables
data.prevPupil      = circshift(data.decision_pupil, 1);
data.prevResp       = circshift(data.resp, 1);
data.prevStim       = circshift(data.stim, 1);
data.resp(data.resp == -1) = 0; % predict response identity

% don't use trials that are at the beginning of each block
trlDif = [0; diff(data.trialnr)];

removeTrls = false(size(trlDif));
removeTrls(trlDif < 1) = true;
removeTrls(trlDif > 1) = true;
removeTrls(find(trlDif > 1) + 1) = true;

%% random effects to include per subject and session
% important, enter all fixed first and then all random effects second in the formula for
% matlab to parse it properly

FullHistMdl = fitglme(data, 'resp ~ 1 + motionstrength + prevResp + (1 + motionstrength + prevResp|newsession)', ...
    'distribution', 'binomial', 'FitMethod', 'Laplace');

% when we allow a different psychometric function (bias + slope) at each
% session per sj, does the fixed effect of history * pupilMod do a better
% job of predicting the data?

% Q1: dummy vs effect code (how to interpret the beta for main effect of
% prevResp when interaction with prevPupil is included?)

FullPupilMdl = fitglme(data, 'resp ~ 1 + motionstrength + prevResp*prevPupil + (1 + motionstrength + prevPupil*prevResp|newsession)', ...
    'distribution', 'binomial', 'FitMethod', 'Laplace')

% important, to do model comparison we need the LaPlace fitting method!
resultsPupil = compare(FullHistMdl, FullPupilMdl, 'CheckNesting', true)
% report as:
% deltaAIC = AIC(mdl1) - AIC(mdl2), chiSq(deltaDF) = LRstat, p = pValue
% see Boehm et al. Neuroimage (response caution paper)

if 0,
    %% split by session to see how robust these results are
    
    for session = 1:3,
        switch session
            case 1
                thisdata = data(data.sessionnr < 3, :);
            case 2
                thisdata = data(data.sessionnr < 5 & data.sessionnr > 2, :);
            case 3
                thisdata = data(data.sessionnr > 5, :);
        end
        
        HistMdl{session} = fitglme(thisdata, 'resp ~ 1 + motionstrength + prevResp + (1 + motionstrength + prevResp|subjnr)', ...
            'distribution', 'binomial', 'FitMethod', 'Laplace');
        
        PupilMdl{session} = fitglme(thisdata, 'resp ~ 1 + motionstrength + prevResp*prevPupil + (1 + motionstrength + prevResp*prevPupil|subjnr)', ...
            'distribution', 'binomial', 'FitMethod', 'Laplace');
        
    end
    
    save('~/Data/UvA_pupil/GrandAverage/glmm_everything.mat');
    
    %%
    load('~/Data/UvA_pupil/GrandAverage/glmm_everything.mat');
    
    clc
    for session = 1:3,
        disp(session)
        tmp = PupilMdl{session}.Coefficients.pValue;
        pvalPup(session) = tmp(end);
        pupilSignific = compare(HistMdl{session}, PupilMdl{session}, 'CheckNesting', true)
    end
    
    % session 3 has a really strong effects but the others not so much...
    
end


