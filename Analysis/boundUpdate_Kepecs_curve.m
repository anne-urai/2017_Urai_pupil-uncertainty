function boundUpdate_Kepecs_curve
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
% See email 18/02/2017
%
% Anne Urai, 2017
% anne.urai@gmail.com

clearvars -except mypath; close all;
global mypath;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

% warnings
warning('error', 'stats:glmfit:PerfectSeparation');
warning('error', 'stats:glmfit:IterationLimit');
warning('error', 'stats:glmfit:IllConditioned');

% get all data
colormap(cbrewer('div', 'RdBu', 64));
set(groot, 'defaultaxescolororder', viridis(3), 'DefaultAxesFontSize', 8);


subjects = 1:27;
for sj = subjects,
    
    % get data
    data     = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % make ME evenly distributed over trials
    for s = unique(data.sessionnr)',
        data.motionstrength(s == data.sessionnr) = ...
            zscore(data.motionstrength(s == data.sessionnr));
    end
    
    % add some stuff
    data.prevcorrect            = circshift(data.correct, 1);
    data.prevresp               = circshift(data.resp, 1);
    data.prevmotionstrength     = circshift(data.motionstrength, 1);
    
    % remove non-continuous trials?
    
    correct = [1 0];
    choices = [-1 1];
    for prevcorrect = 1:2,
        for prevchoice = 1:2,
            
            % subset of trials
            trls = find(data.prevcorrect == correct(prevcorrect) & ...
                data.prevresp == choices(prevchoice));
            thisdat = data(trls, :);
            
            prevdifficulty = findgroups(discretize(abs(thisdat.prevmotionstrength), ...
                [-inf quantile(abs(thisdat.prevmotionstrength), 3) inf]));
            for prevdiff = 1:3,
                thisdat2 = thisdat(prevdifficulty == prevdiff, :);
                
                % fit curve
                b = glmfit(thisdat2.motionstrength, (thisdat2.resp > 0), 'binomial', 'link', 'logit');
                
                % save
                grandavg(sj, prevcorrect, prevchoice, prevdiff, :) = b;
                
            end
        end
    end
end

% =================================================== %
% grand average
% =================================================== %

clf;
grandavg_mean = squeeze(mean(grandavg));
xvals = linspace(-2, 2, 100);
colors = cbrewer('div', 'BrBG', 7); colors = colors([1 2 3 5 6 7], :);
xpos = [3 2 1 4 5 6];

for sp = 1:2, % previous correct vs error
    subplot(3,3,sp); cnt = 1;
    
    for prevchoice = 1:2,
        for prevdiff = 1:3,
            y = glmval(squeeze(grandavg_mean(sp, prevchoice, prevdiff, :)), ...
                xvals, 'logit');
            p(cnt) = plot(xvals, y, 'color', colors(xpos(cnt), :)); hold on;
            cnt = cnt + 1;
        end
    end
    
    switch sp
        case 1
            title('Previous correct');
        case 2
            title('Previous error');
    end
    xlabel('Sensory evidence (z)');
    ylabel('# choice 1');
    % xlim([-0.5 0.5]);
    axis tight; ylim([0 1]);
    set(gca, 'ytick', [0 0.5 1]);
    box off; offsetAxes;
    
end
labtxt = {'choice -1, hard', 'choice -1, medium', 'choice -1, easy', ...
    'choice 1, hard', 'choice 1, medium', 'choice 1, easy'};
l = legend(p(xpos), labtxt(xpos));
l.Position(1) = l.Position(1) + 0.25;
l.Box = 'off';
l.FontSize = 7;
text(3.5, 0.95, 'Previous trial');

% =================================================== %
% only the bias term
% =================================================== %

for sp = 1:2, % previous correct vs error
    subplot(3,3,sp+3); hold on; cnt = 1;
    plot(1:length(labtxt), zeros(size(labtxt)), 'color', [0.5 0.5 0.5], 'linewidth', 0.1);
    for prevchoice = 1:2,
        for prevdiff = 1:3,
            
            p = ploterr(xpos(cnt), nanmean(grandavg(:, sp, prevchoice, prevdiff, 1)), ...
                [], nanstd(grandavg(:, sp, prevchoice, prevdiff, 1)) ./ sqrt(27), 'o', 'abshhxy', 0);
            set(p(1), 'markerfacecolor', colors(xpos(cnt), :), 'markeredgecolor', 'w');
            set(p(2), 'color', colors(xpos(cnt), :));
            
            cnt = cnt + 1;
        end
    end

    % xlabel('Sensory evidence (z)');
    ylabel('Current choice bias');
    % xlim([-0.5 0.5]);
    box off; axis tight;
    ylim([-0.3 0.3]); set(gca, 'ytick', [-0.3 0 0.3]);
    offsetAxes;
    set(gca, 'xtick', 1:cnt-1, 'xticklabel', [-3 -2 -1 1 2 3]);
    xlabel('Previous trial evidence');
end

print(gcf, '-dpdf', sprintf('%s/Figures/boundUpdate_curve.pdf', mypath));

end
