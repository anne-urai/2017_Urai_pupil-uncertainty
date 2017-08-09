function [] = fruendKernels(whichmodulator, field)
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
% determine the subjects based on their plain weights
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator));

if strcmp(whichmodulator, 'plain') && strcmp(field, 'response'),
    % save colors
    colors = cbrewer('div', 'PuOr', 256);
    %colors(1:128, :) = flipud(colors(129:end, :));
    % flip the halves
    colors(128-40:128+40, :) = []; % remove white in the middle
    colors = colors + 0.1; % make a bit paler
    colors(colors > 1) = 1;
    
    for sj = 1:27,
        % find which color this is
        colspace =  linspace(-0.5, 0.5, size(colors, 1));
        colidx      = dsearchn(colspace', dat.(field)(sj, 1));
        mycolmap(sj, :) = colors(colidx, :);
    end
    save(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath), 'mycolmap');
elseif strcmp(field, 'response_pupil'),
    mycolmap = 0.8*ones(27, 3);
else
    % get the colors from all the data combined
    load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));
end

hold on;
for sj = 1:27,
    plot(dat.(field)(sj, :)', 'color', mycolmap(sj, :), 'linewidth', 0.5);
end
s = scatter(ones(1, 27), dat.(field)(:, 1), 10, mycolmap, 'filled', 'linewidth', 0.1, 'markeredgecolor', 'w');

if strcmp(whichmodulator, 'plain') ,
    
    % add the group
    %     [ax, p1, p2] = plotyy(1:7, nanmean(dat.(field)),  ...
    %         1:7, mean(abs(dat.(field))));
    p1 = plot(1:7, nanmean(dat.(field)));
    
    set(gca, 'ycolor', 'k', 'xtick', 1:7, 'ytick', [-0.4 0 0.4], 'ylim', [-.45 .4], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
    % set(ax(2), 'ycolor', [0.4 0.4 0.4], 'xtick', 1:7, 'ytick', [0 0.2], 'ylim', [0 .4], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
    
    set(p1, 'color', 'k', 'linewidth', 1);
    % set(p2, 'color', [0.4 0.4 0.4], 'linewidth', 1);
    axis square; %axis(ax(2), 'square');
    
    if strcmp(field, 'response'),
        
        ylabel(gca,'Choice weight') % label left y-axis
      %  ylabel(ax(2),'|Choice weight|') % label right y-axis
        
    elseif strcmp(field, 'stimulus'),
          ylabel(gca,'Stimulus weight') % label left y-axis
     %   ylabel(ax(2),'|Stimulus weight|') % label right y-axis
        
    end

elseif strcmp(field, 'response_pupil'),
    
    % indicate significance for each lag
    clear h;
    for s = 1:size(dat.(field), 2),
        h(s) = ttest(dat.(field)(:, s));
    end
    h(h < 1) = NaN;
    plot(1:7, -0.2*h, 'k.', 'markersize', 5);
    
end
xlabel('Lags');
offsetAxes;

end
