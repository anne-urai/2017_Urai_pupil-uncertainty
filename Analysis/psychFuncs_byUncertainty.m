function [b] = PsychFuncs_byUncertainty(field)
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
nbins = 6;

% can try this also with all subjects
alldata = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
subjects = unique(alldata.subjnr)';

grandavg.ev     = nan(length(subjects), 2, nbins);
grandavg.acc    = nan(length(subjects), 2, nbins);

for sj = subjects,

    data = alldata(alldata.subjnr == sj, :);

    % project RT out of the pupil and vice versa
    switch field
        case 'rt'
            data.(field) = projectout(zscore(data.rtNorm), zscore(data.decision_pupil));
        case 'decision_pupil'
            data.(field) = projectout(data.(field), zscore(data.rtNorm));
            case 'decision_pupil_only'
            data.(field) = data.decision_pupil;
    end

    % split by low and high RT
    rtMed = quantile(data.(field), 2);

    for r = 1:2,
        switch r
            case 1
                trls = find(data.(field) < rtMed(1));
            case 2
                trls = find(data.(field) > rtMed(2));
        end
        [grandavg.ev(sj, r, :), grandavg.acc(sj, r, :)] = ...
            divideintobins(abs(data.motionstrength(trls)), data.correct(trls), nbins);

    end

    % split by low and high RT
    rtMed = quantile(data.(field), 2);

    for r = 1:2,
        switch r
            case 1
                trls = find(data.(field) < rtMed(1));
            case 2
                trls = find(data.(field) > rtMed(2));
        end

        % fit cumulative weibull to this psychfunc
        data.intensity = data.motionstrength;
        [slope, threshold, lapse] = fitWeibull(abs(data.motionstrength(trls)), data.correct(trls));
        grandavg.b(sj, r, :) = [slope threshold lapse];

    end
end

colors(1,:) = [0.5 0.5 0.5];
colors(2,:) = [0.2 0.2 0.2];

% PLOT
hold on;
for r = 1:2,

    % show the mean curve for the Weibull fit
    newx = linspace(0.1, 5.5, 100);
    plot(newx, 100* weibull(squeeze(nanmean(grandavg.b(:, r, :))), newx), 'color', colors(r, :)) ;

    % datapoints on top
    h = ploterr( squeeze(nanmean(grandavg.ev(:, r, :))), squeeze(100* nanmean(grandavg.acc(:, r, :))), ...
        squeeze( nanstd(grandavg.ev(:, r, :))) ./ sqrt(length(subjects)), ...
        100 * squeeze(nanstd(grandavg.acc(:, r, :))) ./ sqrt(length(subjects)), ...
        '.', 'abshhxy', 0);
    set(h(1), 'color', colors(r,:), 'marker', '.', 'markersize', 12);
    set(h(2), 'color', colors(r,:));
    set(h(3), 'color', colors(r,:));
    handles{r} = h(1);
end

set(gca, 'xlim', [-0.5 5.5], 'ylim', [45 100], 'ytick', 50:25:100);
xlim([-0.5 5.5]); set(gca, 'xtick', [0 2.75 5.5], 'xticklabel', {'weak', 'medium', 'strong'}, 'xminortick', 'off');
ylabel('Accuracy (%)'); xlabel('Evidence');
axis square; box off;

switch field
    case 'rt'
        l = legend([handles{:}], {'fast', 'slow'}, 'location', 'southeast');
    case 'decision_pupil'
        l = legend([handles{:}], {'low', 'high'}, 'location', 'southeast');
    otherwise
        l = legend([handles{:}], {'low', 'high'}, 'location', 'southeast');
end
legend boxoff;
lpos = get(l, 'position');
lpos(1) = lpos(1) + 0.05;
set(l, 'position', lpos);

% take only the second parameter from the Weibull fit (threshold)
b = grandavg.b(:, :, 2);

end
