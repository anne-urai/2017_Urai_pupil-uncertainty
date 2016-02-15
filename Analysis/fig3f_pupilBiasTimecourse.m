function [] = fig3f_pupilBiasTimecourse()
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

clearvars -except mypath; subjects = 1:27;
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens

% get all data
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));

% append all the mean timecourses per condition
for sj = unique(subjects),
    thistable = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % what is this subject's overall bias?
    y = thistable.resp; y(y==-1) = 0;
    [beta,dev,stats] = glmfit(thistable.motionstrength, [y ones(size(y))], ...
        'binomial','link','logit');
    pse = beta(1); % their point of subjective equality
    
    % find their individual 70% correct threshold
    xrange          = abs(thistable.motionstrength);
    yrange          = thistable.correct;
    [newx, newy]    = divideintobins(xrange, yrange, 10);
    % cumulative normal at 0.5
    ft              = fittype( @(b, s, x) 0.5 * (1 + erf( (abs(x-b)/(s * sqrt(2))))));
    fitobj          = fit(xrange, yrange, ft);
    % plot(fitobj); hold on; plot(newx, newy, 'o');
    % fit this to these data
    sortedx         = sort(xrange);
    yval            = feval(fitobj, sortedx);
    threshold       = sortedx(dsearchn(yval, 0.6));
    
    resps = [-1 1];
    cnt  = 0;
    for r = 1:2,
        
        % find trials that were either in line or against the bias
        % first one = same side of pse, second one = other side
        % take only difficult trials
        trls = find(thistable.resp == resps(r) * sign(pse) & abs(thistable.motionstrength) < threshold & thistable.correct == 0);
        pupilchan = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'EyePupil')==1);
        
        cnt = cnt + 1;
        % get all timelock
        alltimelock(sj, cnt, :) = cat(2, squeeze(nanmean(pupilgrandavg.timelock{sj}(1).lock.trial(trls, pupilchan, :)))', ...
            squeeze(nanmean(pupilgrandavg.timelock{sj}(2).lock.trial(trls, pupilchan, :)))', ...
            squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(trls, pupilchan, :)))', ...
            squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(trls, pupilchan, :)))');
    end
end

% color scheme
cols = cbrewer('qual', 'Set1', 4);
cols = cols([2 4], :);

% plot
cla;
ph = boundedline(1:size(alltimelock, 3), squeeze(nanmean(alltimelock)), ...
    permute(squeeze(nanstd(alltimelock)) / sqrt(length(subjects)), [2 3 1]), 'cmap', cols);
ylabel({'Pupil response (z)'});
axis tight;
ph2 = plot(1:2, mean(get(gca, 'ylim'))*ones(2, 10), '.w');
lh = legend(ph2); % make t
lh.String = {'\color[rgb]{0.215686274509804,0.494117647058824,0.721568627450980} with bias', ...
    '\color[rgb]{0.596078431372549,0.305882352941177,0.639215686274510} against bias'};
lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .15;
set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 6);

axis tight; set(gca, 'ytick', [0:0.5:1], 'ylim', [-0.2 1]);
% subfunction to put lines and xlabels at the right spots
plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
    pupilgrandavg.timelock{1}(2).lock, [0], ...
    pupilgrandavg.timelock{1}(3).lock, [0 1], ...
    pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
xlabel('Time (ms)');

%% do stats

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

nsubj = 27;
design = zeros(2,2*nsubj);
for i = 1:nsubj,  design(1,i) = i;        end
for i = 1:nsubj,  design(1,nsubj+i) = i;  end
design(2,1:nsubj)         = 1;
design(2,nsubj+1:2*nsubj) = 2;

cfgstats.design   = design;
cfgstats.uvar     = 1;
cfgstats.ivar     = 2;

for r = 1:2,
    data(r).time       = 1:size(alltimelock, 3);
    data(r).label      = {'EyePupil'};
    data(r).dimord     = 'subj_chan_time';
    data(r).individual = alltimelock(:, r, :);
end

stat      = ft_timelockstatistics(cfgstats, data(1), data(2));
yval      = -0.10;
plot(find(stat.mask==1), yval*ones(1, length(find(stat.mask==1))), '.', 'color', [0.5 0.5 0.5], 'markersize', 4);

end
