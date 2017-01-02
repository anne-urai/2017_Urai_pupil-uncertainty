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

clearvars -except mypath;
global mypath; close all; clc;
mods = {'pupil', 'rt'};
nbins = 3;
for m = 1:2,
    subplot(4,4,m);
    grandavg{m} = postPupilBehaviour(mods{m}, 3, []);
    plot([1 nbins], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;
    
    for r = 1:2,
        x = [1:3] + (r-1)*0.1;
        y = squeeze(grandavg{m}.logisticHistory(:, r, :));
        colors = cbrewer('qual', 'Set2', 3);
        switch r
            case 1
                thismarker = '.';
                thismarkersize = 14;
                thiscolor = colors(1,:);
            case 2
                thismarker = '.';
                thismarkersize = 14;
                thiscolor = colors(3,:);
        end
        
        h = ploterr(x, squeeze(nanmean(y)), [], ...
            squeeze(nanstd(y)) ./sqrt(length(y)), 'k-',  'abshhxy', 0);
        set(h(1), 'color', thiscolor, 'markersize', thismarkersize, 'marker',thismarker);
        set(h(2), 'color', thiscolor); % line color
    end
    
    switch mods{m}
        case 'pupil'
            xticklabs       = repmat({' '}, 1, nbins);
            xticklabs{1}    = 'low';
            xticklabs{end}  = 'high';
            if nbins == 3, xticklabs{2} = 'med'; end
        case 'rt'
            xticklabs       = repmat({' '}, 1, nbins);
            xticklabs{1}    = 'fast';
            xticklabs{end}  = 'slow';
            if nbins == 3, xticklabs{2} = 'med'; end
    end
    
    set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
        'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
    axis square; xlim([0.5 nbins+0.5]);
    ylabel('P(choice = 1');
    
    % do statistics, repeated measures anova across bins
    switch mods{m}
        case 'pupil'
            xlabel('Previous trial pupil');
            set(gca, 'ylim', [0.45 0.55], 'ytick', [0.45:0.05:0.55]);
        case 'rt'
            xlabel('Previous trial RT');
            set(gca, 'ylim', [0.44 0.56], 'ytick', [0.44:0.05:0.56]);
        otherwise
            xlabel(sprintf('Previous trial %s', mods{m}));
    end
    
end

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS5.pdf', mypath));
