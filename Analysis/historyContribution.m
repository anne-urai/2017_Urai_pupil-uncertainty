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

%% variance explained by history as a function of stimulus strength
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plainCoh'));
% dat.variance.stimuli    = reshape(dat.variance.stimuli, [27 5 5]);
% dat.variance.explained  = reshape(dat.variance.explained, [27 5 5]);

% first, average across sessions for each person
% http://stackoverflow.com/questions/36965637/matlab-compute-average-value-corresponding-to-each-column-number-in-matrix
dat.newvariance.stimuli     = [0.625 1.25 2.5 5 10 20 30];
dat.newvariance.explained   = nan(27, length(dat.newvariance.stimuli));
for sj = 1:27,
    A = [dat.variance.stimuli(sj, :)' dat.variance.explained(sj, :)'];
    [Aunq,~,Aind] = unique(A(:,1));
    B = [Aunq,accumarray(Aind,A(:,2),[],@mean)] * 100; % percentage of motion coherence and of history contribtuion
    dat.newvariance.explained(sj, ismember(dat.newvariance.stimuli, B(:, 1))) = B(:, 2);
end

boundedline(dat.newvariance.stimuli, ...
    squeeze(nanmean(dat.newvariance.explained)), ...
    squeeze(nanstd(dat.newvariance.explained)) ./ sqrt(sum(~isnan(dat.newvariance.explained))), ...
    'cmap', [0 0 0]);

axis square; box off;
ylim([-8 100]); xlim([-2 30]);
xlabel('Evidence strength'); ylabel('History contribution (%)');
