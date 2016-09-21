function PsychometricFunction(session)
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
subjects = 1:27;

grandavg.accuracy = nan(length(subjects), nbins);
grandavg.rt       = nan(length(subjects), nbins);

for sj = subjects,
    
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % for supplement
    if isnumeric(session),
        data = data(data.sessionnr == session, :);
    end

    % evidence strength
    data.motionstrength = (abs(data.motionstrength));
    

    % bin - this will return MEAN, not MEDIAN RTs
    [grandavg.ev(sj, :), grandavg.accuracy(sj, :)] = ...
        divideintobins(abs(data.motionstrength), data.correct, nbins, [], @nanmean);
    [~, grandavg.rt(sj, :)] = ...
        divideintobins(data.motionstrength, data.rt, nbins, [], @nanmedian);

end

% plot psychometric func and chronometric func in 1
[ax, hLine1, hLine2] = plotyy(nanmean(grandavg.ev), 100*squeeze(nanmean(grandavg.accuracy)), ...
    1:nbins, squeeze(nanmean(grandavg.rt)));
hLine1.LineStyle = 'none';
hLine2.LineStyle = 'none';

cols(1, :) = [0.2 0.2 0.2];
cols(2, :) = [0.6 0.6 0.6];

% add errorbars
hold(ax(1), 'on');
errorbar(ax(1), nanmean(grandavg.ev), 100*squeeze(nanmean(grandavg.accuracy)), ...
    100*squeeze(nanstd(grandavg.accuracy)) ./ sqrt(length(subjects)), ...
    'color', cols(1,:), 'linewidth', 1);

set(ax(1), 'xlim', [0 5.5], 'xtick', [0.5 3 5.5], 'ylim', [45 101], ...
    'ytick', 50:25:100, 'box', 'off', 'ycolor', ...
    cols(1,:), 'xticklabel', {'weak', 'medium', 'strong'}); 
xlabel(ax(1), 'Evidence');
ylabel(ax(1), 'Accuracy (%)');
axis(ax(1), 'square');

hold(ax(2), 'on');
errorbar(ax(2), nanmean(grandavg.ev), squeeze(nanmean(grandavg.rt)), ...
    squeeze(nanstd(grandavg.rt)) ./ sqrt(length(subjects)), 'color', cols(2,:), 'linewidth', 1);

set(ax(2), 'xlim', [0 5.5], 'xtick', [0.5 3 5.5], 'ylim', [0.23 0.56], ...
    'ytick', 0.25:0.15:0.55, 'box', 'off', 'ycolor', cols(2,:)); 

xlabel(ax(2), 'Evidence');
ylabel(ax(2), 'Reaction time (s)');
axis(ax(2), 'square');

ax(2).YLabel.Rotation = 270;
ypos = ax(2).YLabel.Position;
ypos(1) = ypos(1)+0.8;
ax(2).YLabel.Position = ypos;

axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
  axes(a).FontSize = 7;
end

end