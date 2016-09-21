function [] = decisionStrategies(whichmodulator, errorbars)
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

if ~exist('errorbars', 'var'); errorbars = 1; end
if ~exist('showexample', 'var'); showexample = 1; end

global mypath;

hold on;
plot([-0.8 0.8], [-0.8 0.8], 'color', 'k', 'linewidth', 0.5);
plot([-0.8 0.8], [0.8 -0.8], 'color', 'k', 'linewidth', 0.5);

load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator));

if errorbars,
    for sj = 1:27,
        h = ploterr(dat.response(sj, 1), dat.stimulus(sj, 1), ...
            {dat.responseCI(sj, 1, 1) dat.responseCI(sj, 1, 2)}, ...
            {dat.stimulusCI(sj, 1, 1) dat.stimulusCI(sj, 1, 2)}, '.', 'abshhxy', 0);
        set(h(1), 'color', mycolmap(sj, :), 'markerfacecolor', mycolmap(sj, :));
        set(h(2), 'color', mycolmap(sj, :), 'linewidth', 0.5);
        set(h(3), 'color', mycolmap(sj, :), 'linewidth', 0.5);
    end
else
    for sj = 1:27,
        h = ploterr(dat.response(sj, 1), dat.stimulus(sj, 1), ...
            [], ...
            [], '.', 'abshhxy', 0);
        set(h(1), 'color', mycolmap(sj, :), 'markerfacecolor', mycolmap(sj, :), 'markersize', 10);
    end
end

fz = 6;
text(0, 0.36, 'win stay', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, 0.32, 'lose switch', 'horizontalalignment', 'center', 'fontsize', fz);

text(0, -0.32, 'win switch', 'horizontalalignment', 'center', 'fontsize', fz);
text(0, -0.36, 'lose stay', 'horizontalalignment', 'center', 'fontsize', fz);

text(0.36, .05, 'stay', 'rotation', 270, 'fontsize', fz);
text(-0.36, -.06, 'switch', 'rotation', 90, 'fontsize', fz);

% layout
maxlim = 0.5;
xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
maxlim = 0.4;
set(gca, 'xtick', -maxlim:maxlim:maxlim, 'ytick', -maxlim:maxlim:maxlim);
xlabel('Choice weight'); ylabel('Stimulus weight');
box on; axis square;

end