function [] = motionEnergy_Probability()
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

if 0,
    % for each sj, check if this worked
    for sj = 1:27,
        t = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        subplot(5,6,sj);
        plot(t.stim .* t.coherence, t.motionstrength, '.');
        title(sprintf('P%02d', sj)); axis tight;
        box off;
        ylim([-7 7]);
    end
end

%% group level results
t = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
stim = t.coherence;
stim = round((stim) * 100000)/100000;
stim((abs(stim-0.01)) < 0.000001) = 0.0125;
stim = t.stim .* stim;
strength = t.motionstrength;

% !!! add zero stimlevels into image MAT to ensure linearity
stepsize = min(diff(unique(stim)));
stimlevels = min(unique(stim)):stepsize:max(unique(stim));
stimlevels = [-0.31 stimlevels 0.31];
n = zeros(length(stimlevels), 99);

% find
for s = 1:length(stimlevels),
    trls = find(abs(stim - stimlevels(s)) < 0.000001);
    edges = linspace(min(strength), max(strength), 100);
    [n(s, :), edges] = histcounts(strength(trls), edges, ...
        'Normalization', 'probability');
end

% check
assert(isequal(length(unique(stim)), length(find(nansum(n, 2) > 0))), 'spacing not quite right');
colormap(inferno);
imagesc(stimlevels, edges, n');
axis square;
set(gca, 'ydir', 'normal');
xlabel({'Motion coherence (%)'});
ylabel({'Motion energy (a.u.)'});

%offsetAxes
set(gca, 'tickdir', 'out', 'box', 'off');
set(gca, 'xcolor', 'k', 'ycolor', 'k');
set(gca, 'xtick', [-0.3 0 0.3], 'xticklabel', [-30 0 30], 'xticklabelrotation', 0);
set(gca, 'ytick', -6:6:6);
set(gca, 'XAxisLocation','top');

axpos = get(gca, 'Position');
drawnow;

% add the colorbar, make it prettier
handles = colorbar('Location', 'EastOutside');
handles.TickDirection = 'out';
handles.Box = 'off';
%handles.Label.String = 'Probability';

% this looks okay, but the colorbar is very wide. Let's change that!
% get original axes
cpos = handles.Position;
cpos(3) = 0.5*cpos(3);
cpos(1) = cpos(1) + 0.05; % to the right
handles.Position = cpos;
drawnow;

% restore axis pos
set(gca, 'position', axpos);
drawnow;

end

