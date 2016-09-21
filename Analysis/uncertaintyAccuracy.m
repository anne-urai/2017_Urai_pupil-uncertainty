function [b] = uncertaintyAccuracy(field)
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

global mypath;
nbins = 12;

% can try this also with all subjects
subjects = 1:27;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grandavg.rt     = nan(length(subjects), nbins);
grandavg.acc    = nan(length(subjects), nbins);
grandavg.b      = nan(length(subjects), 2);

for sj = subjects,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % bin
    [grandavg.rt(sj, :), grandavg.acc(sj, :)] = divideintobins(data.(field), data.correct, nbins);
    
    % project RT out of the pupil and vice versa
    switch field
        case 'rt'
            data.(field) = projectout(zscore(data.rtNorm), zscore((data.decision_pupil)));
        case 'decision_pupil'
            data.(field) = projectout(zscore(data.decision_pupil), zscore(data.rtNorm));
    end
    
    % also do logistic regression
    grandavg.b(sj, :) = glmfit(zscore(data.(field)), data.correct, ...
        'binomial','link','logit');
end

h = boundedline(1:nbins, 100* nanmean(grandavg.acc), ...
    100 * nanstd(grandavg.acc) ./ sqrt(length(subjects)), 'cmap', [0 0 0]);

% stats
xlim([0 nbins]); set(gca, 'xtick', [1 nbins/2 nbins], 'xminortick', 'off');
switch field
    case 'rt'
        ylim([45 90]); set(gca, 'ytick', [50:20:100]);
        xlabel('Reaction time');
        set(gca,  'xticklabel', {'fast', 'med', 'slow'});
    case 'decision_pupil'
        xlabel('Pupil response');
        ylim([68 80]); set(gca, 'ytick', [70:10:100]);
        set(gca,  'xticklabel', {'low', 'medium', 'high'});
end
ylabel('Accuracy (%)');
axis square;

b = grandavg.b;
end