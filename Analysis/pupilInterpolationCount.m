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

%% count the nr of trials with > 50% of samples interpolated
global mypath;

for sj = 1:27,
    
    load(sprintf('%s/Data/P%02d_alleyeNaN.mat', mypath, sj));
    
    percentageNaN = @(x) mean(isnan(x), 2);
    % for each trial, get the percentage of interpolated samples
    interpolatedSamples = cellfun(percentageNaN, data.trial, 'uniformoutput', 0);
    interpolatedSamples = cat(2, interpolatedSamples{:});
    interpolatedSamples = interpolatedSamples(find(strcmp(data.label, 'EyePupil')==1), :); % take the pupil chan
    
    % for each subject, compute the percentage of trials with > 50% rejected samples
    grandavg.rejectTrials(sj) = 100 * mean(interpolatedSamples > 0.5);
    
    % for computing the total
    grandavg.alltrials(sj) = length(data.trial);
    grandavg.nrreject(sj) = sum(interpolatedSamples > 0.5);
end

% total
total = 100 * sum(grandavg.nrreject) ./ sum(grandavg.alltrials);
fprintf('\n\nrange of trials with > 50%% interpolation: %.3f %% to %.3f %%, mean %.3f %%, total %.3f %% \n', ...
    min(grandavg.rejectTrials), max(grandavg.rejectTrials), mean(grandavg.rejectTrials), total);

