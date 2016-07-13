function data = addUncertainty2data

%% another test, correlation matrix
close all; clear; clc;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

% for sj = 1:27,
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));

% fit a probit slope so we can get the sigma
b       = glmfit(data.motionstrength, (data.resp > 0), 'binomial', 'link', 'probit');
sigma   = 1/b(2);  % standard deviation at these values is the inverse!
bound   = -b(1); % the bound is negative if people say

% for each trial, compute the average level of uncertainty
data.uncertainty    = arrayfun(@simulateUncertainty, abs(data.motionstrength), ...
    data.correct, sigma*ones(length(data.correct), 1), bound*ones(length(data.correct), 1));

% absolute evidence strength
data.evidence       = abs(data.motionstrength);

% add a measure of repeat
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

% remove the repeat measures for trials with non-consecutive trial nrs
trlDif                      = [diff(data.trialnr); 0];
removeTrls                  = false(size(trlDif));
removeTrls(trlDif < 1)      = true;
removeTrls(trlDif > 1)      = true;
data.stimrepeat(removeTrls) = NaN;
data.repeat(removeTrls)     = NaN;
end