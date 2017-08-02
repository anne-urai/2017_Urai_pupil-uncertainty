function [] = pupilUncertaintyTimecourse(plotAll)
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
if ~exist('plotall', 'var'); plotAll = false; end % in principle, only the main evidence strength regressor

subjects = 1:27;
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens
warning('error', 'stats:regress:RankDefDesignMat'); % stop if this happens

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% instead of binning, use regression beta across samples
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = cbrewer('qual', 'Set1', 9);
colors = colors([1 2 9], :); % red blue grey

if exist(sprintf('%s/Data/GrandAverage/pupilRegressionBetas.mat', mypath), 'file'), ...
        load(sprintf('%s/Data/GrandAverage/pupilRegressionBetas.mat', mypath)); % dont redo, takes time;
else
    for sj = unique(subjects),
        
        pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(3).lock.label, 'EyePupil')==1);
        
        % get all timelock
        tl = cat(2, ...
            squeeze(pupilgrandavg.timelock{sj}(1).lock.trial(:, pupilchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(2).lock.trial(:, pupilchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :)));
        
        % do the same for the gazepos, x and y
        xchan    = find(strcmp(pupilgrandavg.timelock{sj}(3).lock.label, 'EyeH')==1);
        gazex = cat(2, ...
            squeeze(pupilgrandavg.timelock{sj}(1).lock.trial(:, xchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(2).lock.trial(:, xchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(3).lock.trial(:, xchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(4).lock.trial(:, xchan, :)));
        gazex = gazex - nanmean(gazex(:)); % normalize
        
        ychan    = find(strcmp(pupilgrandavg.timelock{sj}(3).lock.label, 'EyeV')==1);
        gazey = cat(2, ...
            squeeze(pupilgrandavg.timelock{sj}(1).lock.trial(:, ychan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(2).lock.trial(:, ychan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(3).lock.trial(:, ychan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(4).lock.trial(:, ychan, :)));
        gazey = gazey - nanmean(gazey(:)); % normalize
        
        % run a separate glm for correct and error
        cors = [0 1];
        signific = ones(2, size(tl, 2));
        
        for c = 1:2,
            
            thistabledat = pupilgrandavg.timelock{sj}(4).lock.trialinfo;
            trls = find(thistabledat(:, 8) == cors(c));
            
            % add zscored log(RT) as another predictor
            designM = [ones(length(trls), 1) zscore(abs(thistabledat(trls, 4))) ...
                zscore(thistabledat(trls, 15))];
            
            % regress for each sample
            for s = 1:size(tl, 2),
                
                % add the x and y position of this sample?
                designM2 = [designM gazex(trls, s) gazey(trls, s)];
                
                % when zscoring each sample's timecourse, the nice intercept shape
                % dissappears (duh)
                [b, bint, ~, ~, stats] = regress((tl(trls, s)), designM2);
                grandavg.beta(find(sj==subjects), c, s, :)    = b;
                grandavg.bint(find(sj==subjects), c, s, :)    = [b - bint(:, 1)];
                grandavg.rsq(find(sj==subjects), c, s, :)     = stats(1);
                
                signific(c, s) = stats(3);
            end
        end
        
    end
    % save for next time
    savefast(sprintf('%s/Data/GrandAverage/pupilRegressionBetas.mat', mypath), 'grandavg', 'signific'); % dont redo, takes time;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the timecourse of regression coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotAll,
    whichBetas2plot = 1:size(designM2, 2);
else
    whichBetas2plot = 2; % only for main figure
end

for whichbeta = whichBetas2plot,
    % line to indicate zero
    plot([0 size(grandavg.beta, 3)], [0 0], '-', 'color', 'k', 'LineWidth', 0.2);
    
    hold on;
    boundedline(1:size(grandavg.beta, 3), squeeze(nanmean(grandavg.beta(:, :, :, whichbeta)))', ...
        permute(squeeze(nanstd(grandavg.beta(:, :, :, whichbeta))) / sqrt(length(subjects)), [2 3 1]), ...
        'cmap', colors);
    axis tight;
    ylim([-0.1 0.1]);
    thisax = gca;
    
    hold on;
    % subfunction to put lines and xlabels at the right spots
    plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
        pupilgrandavg.timelock{1}(2).lock, [0], ...
        pupilgrandavg.timelock{1}(3).lock, [0 1], ...
        pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
    
    % indicate the grey area we use for getting single-trial scalars
    xticks = get(gca, 'xtick');
    finalx = xticks(5);
    startx = xticks(5) - (0.25* (xticks(6)-xticks(5)));
    a = area(startx:finalx-1, ...
        ones(1, finalx-startx) * max(get(gca, 'ylim')), ...
        min(get(gca, 'ylim')));
    a.FaceColor = [0.9 0.9 0.9];
    a.EdgeColor = 'none';
    a.ShowBaseLine = 'off';
    
    if exist(sprintf('%s/Data/GrandAverage/pupilRegressionSignificantCluster.mat', mypath), 'file'), % skip stats for the intercept
        load(sprintf('%s/Data/GrandAverage/pupilRegressionSignificantCluster.mat', mypath));
    else
        % on this, do cluster based permutation testing over the time
        % dimension - because the samples are not temporally independent...
        
        for c = 1:2,
            % put into fieldtrip-style
            grandavg_thiscorr(c).time       = 1:size(grandavg.beta, 3);
            grandavg_thiscorr(c).label      = {'EyePupil'};
            grandavg_thiscorr(c).dimord     = 'subj_chan_time';
            grandavg_thiscorr(c).individual = permute(squeeze(grandavg.beta(:, c, :, whichbeta)), [1 3 2]);
            
            % compare against null
            grandavg_zero = grandavg_thiscorr(c);
            grandavg_zero.individual = zeros(size(grandavg_thiscorr(c).individual));
            
            stat{c} = clusterStat(grandavg_thiscorr(c), grandavg_zero, length(subjects));
        end
        
        % also test their difference
        stat{3} = clusterStat(grandavg_thiscorr(1), grandavg_thiscorr(2), length(subjects));
        savefast(sprintf('%s/Data/GrandAverage/pupilRegressionSignificantCluster.mat', mypath), 'stat');
    end
    
    yval = -0.12;
    % plot
    for c = 1:2,
        % plot on top
        p(c) = plot(find(stat{c}.mask==1), yval*ones(1, length(find(stat{c}.mask==1))), '.', 'color', colors(c, :), 'markersize', 4);
        yval = yval - 0.08*range(get(thisax, 'ylim'));
        
        % output when this starts to become significant
        signific  = find(stat{c}.mask==1);
        alltiming = [pupilgrandavg.timelock{1}(3).lock.time pupilgrandavg.timelock{1}(4).lock.time];
        try
            disp([alltiming(signific(1)) alltiming(signific(end))]);
        end
    end
    p(3) = plot(find(stat{3}.mask==1), yval*ones(1, length(find(stat{3}.mask==1))), '.', 'color', colors(3, :), 'markersize', 4);
    
    
    ylabel('Beta weight (a.u.)');
    ylim([-0.18 0.1]);
    
    hold on;
    ph2 = plot(180:182, mean(get(thisax, 'ylim'))*ones(3, 10), '.w');
    lh = legend(ph2);
    legnames = {'error', 'correct', 'error v correct'};
    for i = 1:length(legnames),
        str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
    end
    lh.String = str;
    lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .1;
    lpos(2) = lpos(2) - 0.08;
    set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 6);
end

end
