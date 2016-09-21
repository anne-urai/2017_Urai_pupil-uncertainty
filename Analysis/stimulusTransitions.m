% This code reproduces the analyses in the paper
% Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven 
% by decision uncertainty and alters serial choice bias. 
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai, 2016
% anne.urai@gmail.com

%% sequential biases in stimuli
stimRep.transition = nan(27, 1);
for sj  = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    stim = data.stim;
    stim(logical([(diff(data.trialnr) ~= 1); 1])) = NaN; % ignore stimuli that cross block boundaries
    stim(stim == -1) = 0;
    stimRep.transition(sj, :) = nanmean(abs(diff(stim)));
end
fprintf('transitionprob = %.3f, range %.3f-%.3f \n', ...
    nanmean(stimRep.transition(:)), min(stimRep.transition(:)), max(stimRep.transition(:)));

% correlate with individual weights
load(sprintf('%s/Data/GrandAverage/historyweights_plain.mat', mypath));
[rho, pval] = corr(stimRep.transition, dat.response(:, 1));
bf10 = corrbf(rho, 27);
fprintf('r = %.3f, p = %.3f, BF10 = %.3f \n', rho, pval, bf10);

% and stimulus weights
[rho, pval] = corr(stimRep.transition, dat.stimulus(:, 1));
bf10 = corrbf(rho, 27);
fprintf('r = %.3f, p = %.3f, BF10 = %.3f', rho, pval, bf10);


%% compute autocorrelation in evidence strength
stimRep.rho = nan(27, 1);
stimRep.pval = nan(27, 1);
for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    stim1       = abs(data.motionstrength);
    stim1(logical([(diff(data.trialnr) ~= 1); 1])) = NaN; % ignore stimuli that cross block boundaries
    [rho, pval] = corr(stim1, circshift(stim1, 1), 'rows', 'complete');
    stimRep.rho(sj) = rho;
    stimRep.pval(sj) = pval;
end
fprintf('rho = %.3f, range %.3f-%.3f, significant in %d out of %d sessions', ...
    nanmean(stimRep.rho(:)), min(stimRep.rho(:)), max(stimRep.rho(:)), ...
    length(find(stimRep.pval(:) < 0.05)), sum(~isnan(stimRep.pval(:))));

% correlate with individual weights
load(sprintf('%s/Data/GrandAverage/historyweights_plain.mat', mypath));
[rho, pval] = corr(stimRep.rho, dat.response(:, 1));
bf10 = corrbf(rho, 27);
fprintf('r = %.3f, p = %.3f, BF10 = %.3f', rho, pval, bf10);

% and stimulus weights
[rho, pval] = corr(stimRep.rho, dat.stimulus(:, 1));
bf10 = corrbf(rho, 27);
fprintf('r = %.3f, p = %.3f, BF10 = %.3f', rho, pval, bf10);

% correlate with individual pupil x choice weights
load(sprintf('%s/Data/GrandAverage/historyweights_pupil+rt.mat', mypath));
[rho, pval] = corr(stimRep.rho, dat.response_pupil(:, 1));
bf10 = corrbf(rho, 27);
fprintf('r = %.3f, p = %.3f, BF10 = %.3f \n', rho, pval, bf10);

% and stimulus weights
[rho, pval] = corr(stimRep.rho, dat.response_rt(:, 1));
bf10 = corrbf(rho, 27);
fprintf('r = %.3f, p = %.3f, BF10 = %.3f', rho, pval, bf10);
