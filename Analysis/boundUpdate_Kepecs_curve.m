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
set(groot, 'defaultaxescolororder', viridis(3), 'DefaultAxesFontSize', 6);
plotIndividual = false;

colors = cbrewer('div', 'BrBG', 9); colors = colors([1:3 end-3:end], :);
colors = cbrewer('div', 'PiYG', 7); colors = colors([1 end], :);

labtxt = {'choice -1', 'choice 1'};

% sort by history weigths
load(sprintf('%s/Data/Grandaverage/historyweights_plain.mat', mypath));
historyweights  = dat.response(:, 1);
[val, spidx]    = sort(historyweights);

for subgroups = 1:3,
    clf;
    switch subgroups
        case 1
            subjects = 1:27;
            groupname = 'allsj';
        case 2
            subjects = find(dat.response(:, 1) < 0)';
            groupname = 'alternators';
        case 3
            subjects = find(dat.response(:, 1) > 0)';
            groupname = 'repeaters';
    end
    spcnt = 1;
    correct = [1 0];
    
    splitgroups = {'coherence', 'pupil', 'RT'};
    for plots = 1:length(splitgroups),
        splitname = splitgroups{plots};
        grandavg = nan(27, 2, 2, 3, 2);
        totalgrandavg = nan(27, 2);
        
        for prevcorrect = 1:2,
            
            for sj = subjects,
                
                % get data
                data     = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
                
                % make ME evenly distributed over trials
                for s = unique(data.sessionnr)',
                    data.motionstrength(s == data.sessionnr) = ...
                        zscore(data.motionstrength(s == data.sessionnr));
                end
                
                % fit curve
                b = glmfit(data.motionstrength, ...
                    (data.resp > 0), 'binomial', 'link', 'logit');
                
                % save
                totalgrandavg(sj, :) = b;
                
                % add some stuff
                data.prevcorrect            = circshift(data.correct, 1);
                data.prevresp               = circshift(data.resp, 1);
                
                switch splitname
                    case 'RT'
                        data.prevmotionstrength     = circshift(zscore(log(data.rt + 0.5)) + 10, 1);
                    case 'coherence'
                        data.prevmotionstrength     = circshift(abs(data.motionstrength), 1);
                    case 'pupil'
                        data.prevmotionstrength     = circshift(data.decision_pupil, 1);
                end
                
                if plotIndividual,
                    % plot
                    subplot(5,6,spidx(sj));
                    hold on; cnt = 1;
                    plot(1:length(labtxt), zeros(size(labtxt)), 'color', [0.5 0.5 0.5], 'linewidth', 0.1);
                end
                
                % remove non-continuous trials?
                choices = [-1 1];
                for prevchoice = 1:2,
                    
                    % subset of trials
                    trls = find(data.prevcorrect == correct(prevcorrect) & ...
                        data.prevresp == choices(prevchoice));
                    thisdat = data(trls, :);
                    
                    % bins of whatever wer are splitting by
                    prevdifficulty = findgroups(discretize(thisdat.prevmotionstrength, ...
                        [-inf quantile(thisdat.prevmotionstrength, 3) inf]));
                    for prevdiff = 1:3,
                        thisdat2 = thisdat(prevdifficulty == prevdiff, :);
                        
                        % fit curve
                        [b, ~, stats] = glmfit(thisdat2.motionstrength, ...
                            (thisdat2.resp > 0), 'binomial', 'link', 'logit');
                        
                        % save
                        grandavg(sj, prevcorrect, prevchoice, prevdiff, :) = b;
                        
                        if plotIndividual
                            % plot
                            p = ploterr(xpos(cnt), b(1), [], stats.se(1), 'o', 'abshhxy', 0);
                            set(p(1), 'markerfacecolor', colors(xpos(cnt), :), 'markeredgecolor', 'w');
                            set(p(2), 'color', colors(xpos(cnt), :));
                            cnt = cnt + 1;
                        end
                        
                    end
                end
                
                if plotIndividual,
                    title(sprintf('P%02d', sj), 'fontweight', 'normal'); box off;
                    ylims = get(gca, 'ylim');
                    ylims = max(abs(ylims));
                    ylim([-ylims ylims]);
                    axis tight; axis square;
                    % offsetAxes;
                    set(gca, 'xtick', 1:cnt-1, 'xticklabel', [-3 -2 -1 1 2 3], ...
                        'xlim', [0.5 max(get(gca, 'xlim'))]);
                end
            end
            
            % =================================================== %
            % normalise by the total bias for every person
            % =================================================== %
            
            grandavg = bsxfun(@minus, grandavg, totalgrandavg(:, 1));
            
            % =================================================== %
            % grand average
            % =================================================== %
            
            subplot(4,4,spcnt); spcnt = spcnt + 1;
            hold on; cnt = 1;
            plot(1:3, zeros(1,3), 'color', [0.5 0.5 0.5], 'linewidth', 0.1);
            for prevchoice = 1:2,
                
                p = ploterr(1:3, squeeze(nanmean(grandavg(:, prevcorrect, prevchoice, :, 1))), ...
                    [], squeeze(nanstd(grandavg(:, prevcorrect, prevchoice, :, 1))) ./ sqrt(length(subjects)), ...
                    'o-', 'abshhxy', 0);
                set(p(1), 'markerfacecolor', colors(cnt, :), 'markeredgecolor', 'w', 'color', colors(cnt, :));
                set(p(2), 'color', colors(cnt, :));
                h(cnt) = p(1);
                cnt = cnt + 1;
            end
            
            box off; axis tight;
            ylim([-0.5 0.5]);
            set(gca, 'ytick', [-0.5 0 0.5]);
            set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'}, ...
              'xlim', [0.1 3.1]);
            
            if plots == 1,
                switch correct(prevcorrect)
                    case 0
                        title('Previous error');
                    case 1
                        title('Previous correct');
                end
            end
            
            xlabel(sprintf('Previous %s bin', splitname));
            ylabel('Current choice bias');
           %  offsetAxes;
            
        end
        spcnt = spcnt + 2;
    end
    % if c == 2,
    l = legend(h, labtxt);
    legend boxoff;
    l.Position(1) = l.Position(1) + 0.2;
    print(gcf, '-dpdf', sprintf('%s/Figures/boundUpdate_total_%s.pdf', mypath, groupname));
end

end
