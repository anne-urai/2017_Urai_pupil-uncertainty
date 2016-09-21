function pupilDprimeCriterionTimecourse(regressOutRT)
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

%{
Pupil bias results, before COSYNE
    - for each sj, get their intrinsic bias (sign)
    - split trials in with-bias and against-bias
    - plot pupil regression timecourses for those, clusterstats on their difference
        - all trials
        - error trials
        - only difficulty trials (threshold <70% correct or so)

    - talk to JW; in window of 250 ms long (move along in steps of 50), bin pupil into 5 bins,
    compute d' and criterion for each of those bins, then correlation
    coefficient between pupil and those two behavioural measures. Then plot this over time!
%}

global mypath;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPLIT into with-bias and against-bias
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except mypath regressOutRT; subjects = 1:27;
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens

% get all data
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));

% append all the mean timecourses per condition
for sj = unique(subjects),
    pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'EyePupil')==1);
    % get all timelock
    alltimelock(sj, 1, :) = cat(2, squeeze(nanmean(pupilgrandavg.timelock{sj}(1).lock.trial(:, pupilchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(2).lock.trial(:, pupilchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :)))');
    
    grandavg(sj).lock = cat(1, squeeze(pupilgrandavg.timelock{sj}(1).lock.trial(:, pupilchan, :))', ...
        squeeze(pupilgrandavg.timelock{sj}(2).lock.trial(:, pupilchan, :))', ...
        squeeze(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, :))', ...
        squeeze(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :))')';
    
    % before doing anything else, take the variance associated with RT out
    % of the pupil signal on each sample!
    if regressOutRT,
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        for i = 1:size(grandavg(sj).lock, 2),
            grandavg(sj).lock(:, i) = projectout(grandavg(sj).lock(:, i), zscore(log(data.rt+0.1)));
        end
    end
    
end

% move in a sliding window of 250 ms, steps of 50
fsample = 100;
timesteps = 13:5:612;
timewinhalf = 12; % so total will be 250 ms, 120 on either side plus the sample itself
nbins = 5;
for sj = subjects,
    
    % get behav data
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data.resp(data.resp == -1) = 0;
    
    for t = timesteps,
        pupilval = mean(grandavg(sj).lock(:, t-timewinhalf:t+timewinhalf), 2);
        
        % divide into 5 bins
        pupQs = quantile(pupilval, nbins-1);
        betas = nan(nbins, 2);
        
        for p = 1:nbins,
            switch p
                case 1
                    trls = find(pupilval < pupQs(p));
                case nbins
                    trls = find(pupilval > pupQs(p-1) );
                otherwise
                    trls = find(pupilval > pupQs(p-1) & pupilval < pupQs(p) );
            end
            
            % based on those trials, get bias and sensitivity
            [beta,dev,stats] = glmfit(data.motionstrength(trls), [data.resp(trls) ones(size(trls))], ...
                'binomial','link','logit');
            betas(p, :) = beta;
            
        end % pupil bins
        
        % correlation between pupil and dprime
        results.bias(sj, find(t==timesteps))        = corr(transpose(1:nbins), abs(betas(:, 1)));
        results.sensitivity(sj, find(t==timesteps)) = corr(transpose(1:nbins), betas(:, 2));
    end
end

effectsize = cat(1, mean(results.bias), mean(results.sensitivity));
effectsem  = cat(1, std(results.bias) ./sqrt(27), std(results.sensitivity)./sqrt(27));

% color scheme
cols = cbrewer('qual', 'Set1', 4);
cols = cols([2 4], :);

% plot
newtimestep = linspace(1, length(squeeze(alltimelock(1, 1, :))), length(timesteps));
plot([newtimestep(1) newtimestep(end)], [0 0],  'k', 'linewidth', 0.2);
hold on;
boundedline(newtimestep, effectsize, ...
    permute(effectsem, [2 3 1]), 'cmap', cols, 'alpha');
axis tight; set(gca, 'ytick', -1:0.2:1);
ph2 = plot(1:2, mean(get(gca, 'ylim'))*ones(2, 10), '.w');
lh = legend(ph2); % make t
lh.String = {'\color[rgb]{0.215686274509804,0.494117647058824,0.721568627450980} bias', ...
    '\color[rgb]{0.596078431372549,0.305882352941177,0.639215686274510} sensitivity', ...
    '\color[rgb]{0.6,0.6,0.6} bias v sensitivity'};
lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .15;
set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 6);

plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
    pupilgrandavg.timelock{1}(2).lock, [0], ...
    pupilgrandavg.timelock{1}(3).lock, [0 1], ...
    pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
ylabel('Pupil correlation (rho)');

% on this, do cluster based permutation testing over the time
% dimension - because the samples are not temporally independent...
yval = -0.95;
thisax = gca;
for c = 1:2,
    % put into fieldtrip-style
    grandavg_thiscorr(c).time       = newtimestep;
    grandavg_thiscorr(c).label      = {'EyePupil'};
    grandavg_thiscorr(c).dimord     = 'subj_chan_time';
    switch c
        case 1
            grandavg_thiscorr(c).individual = permute(squeeze(results.bias), [1 3 2]);
        case 2
            grandavg_thiscorr(c).individual = permute(squeeze(results.sensitivity), [1 3 2]);
    end
    
    % compare against null
    grandavg_zero = grandavg_thiscorr(c);
    grandavg_zero.individual = zeros(size(grandavg_thiscorr(c).individual));
    stat{c} = clusterStat(grandavg_thiscorr(c), grandavg_zero, length(subjects));
    
    % plot on top
    try
        p(c) = plot(newtimestep(find(stat{c}.mask==1)), yval*ones(1, length(find(stat{c}.mask==1))), '.', 'color', cols(c, :), 'markersize', 4);
        yval = yval - 0.08*range(get(thisax, 'ylim'));
    end
end

% also test their difference
stat{3} = clusterStat(grandavg_thiscorr(1), grandavg_thiscorr(2), length(subjects));
p(3) = plot(newtimestep(find(stat{3}.mask==1)), yval*ones(1, length(find(stat{3}.mask==1))), '.', 'color', [0.6 0.6 0.6], 'markersize', 4);
%ylim([-1.2 0.3]);

if regressOutRT,
    title('RT regressed out');
end

end


function stat = clusterStat(data1, data2, nsubj)

% do cluster stats across the group
cfgstats                  = [];
cfgstats.method           = 'montecarlo'; % permutation test
cfgstats.statistic        = 'ft_statfun_depsamplesT'; % dependent samples ttest

% do cluster correction
cfgstats.correctm         = 'cluster';
cfgstats.clusteralpha     = 0.05;
% cfgstats.clusterstatistic = 'maxsize'; % weighted cluster mass needs cfg.wcm_weight...
% cfgstats.minnbchan        = 1; % average over chans
cfgstats.tail             = 0; % two-tailed!
cfgstats.clustertail      = 0; % two-tailed!
cfgstats.alpha            = 0.025;
cfgstats.numrandomization = 1000; % make sure this is large enough
cfgstats.randomseed       = 1; % make the stats reproducible!

% use only our preselected sensors for the time being
cfgstats.channel          = 'EyePupil';

% specifies with which sensors other sensors can form clusters
cfgstats.neighbours       = []; % only cluster over data and time

design = zeros(2,2*nsubj);
for i = 1:nsubj,  design(1,i) = i;        end
for i = 1:nsubj,  design(1,nsubj+i) = i;  end
design(2,1:nsubj)         = 1;
design(2,nsubj+1:2*nsubj) = 2;

cfgstats.design   = design;
cfgstats.uvar     = 1;
cfgstats.ivar     = 2;

stat      = ...
    ft_timelockstatistics(cfgstats, ...
    data1, data2);

end
