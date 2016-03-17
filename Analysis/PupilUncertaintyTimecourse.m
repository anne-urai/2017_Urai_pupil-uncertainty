function [] = PupilUncertaintyTimecourse(plotAll)
% 3. timecourse of regression betas, separately for correct and error

global mypath;
doStats             = true; % permutation statistics across the group - this takes a while!
plotIndividual      = false; % plots all the individual beta timecourses, takes forever
if ~exist('plotall', 'var'); plotAll = false; end % in principle, only the main evidence strength regressor

subjects = 1:27;
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));

warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens
warning('error', 'stats:regress:RankDefDesignMat'); % stop if this happens

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% instead of binning, use regression beta across samples
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = cbrewer('qual', 'Set1', 9);
cols = colors([1 3 9], :); % red green grey

if plotIndividual, figure; end
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
        
        % add log(RT) as another predictor
        designM = [ones(length(trls), 1) zscore(abs(thistabledat(trls, 4))) ...
            zscore(log(thistabledat(trls, 6)+0.1))];
        
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
    
    if plotIndividual,
        % plot this sj
        subplot(6,5,find(sj==subjects));
        boundedline(1:length(signific), squeeze(grandavg.beta(find(sj==subjects), :, :, 2)), ...
            permute(squeeze(grandavg.bint(find(sj==subjects), :, :, 2)), [2 3 1]), 'cmap', cols, 'alpha');
        axis tight;
        set(gca, 'TickDir', 'out', 'XLim', [0 size(grandavg.beta, 3)], 'box', 'off');
        ylims = get(gca, 'ylim');
        plot(find(signific(1, :)<0.05), ylims(1)-0.2*ones(1, length(find(signific(1, :)<0.05))), '.', 'color', cols(1, :));
        plot(find(signific(2, :)<0.05), ylims(1)-0.3*ones(1, length(find(signific(2, :)<0.05))), '.', 'color', cols(2, :));
        
        plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
            pupilgrandavg.timelock{1}(2).lock, [0], ...
            pupilgrandavg.timelock{1}(3).lock, [0 1], ...
            pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
        set(gca, 'xticklabel', []);
    end
end

% save to disk to later do the correlation with learning
save(sprintf('%s/Data/GrandAverage/pupilRegressionBetas.mat', mypath), 'grandavg');

if plotIndividual,
    % subfunction to put lines and xlabels at the right spots
    suplabel('Time (s)', 'x');
    suplabel('Beta (a.u.)', 'y');
    
    if RTstratification,
        suplabel('Beta of coherence with RT as additional predictor', 't');
    else
        suplabel('Beta of coherence', 't');
    end
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
        'cmap', cols);
    axis tight; 
    ylim([-0.1 0.1]); yval = -0.11;
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
    
    if doStats, % skip stats for the intercept
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
            
            % plot on top
            try
                p(c) = plot(find(stat{c}.mask==1), yval*ones(1, length(find(stat{c}.mask==1))), '.', 'color', cols(c, :), 'markersize', 4);
            end
            yval = yval - 0.08*range(get(thisax, 'ylim'));
            
            % output when this starts to become significant
            signific  = find(stat{c}.mask==1);
            alltiming = [pupilgrandavg.timelock{1}(3).lock.time pupilgrandavg.timelock{1}(4).lock.time];
            try
                disp([alltiming(signific(1)) alltiming(signific(end))]);
            end
        end
        
        try
            % also test their difference
            stat{3} = clusterStat(grandavg_thiscorr(1), grandavg_thiscorr(2), length(subjects));
            p(3) = plot(find(stat{3}.mask==1), yval*ones(1, length(find(stat{3}.mask==1))), '.', 'color', cols(3, :), 'markersize', 4);
        end
    end
    
    ylabel('Beta evidence');
    ylim([-0.18 0.1]);
    
    if whichbeta == whichBetas2plot(end) && doStats,
        hold on;
        ph2 = plot(180:182, mean(get(thisax, 'ylim'))*ones(3, 10), '.w');
        lh = legend(ph2); % make t
        lh.String = {'\color[rgb]{0.915294117647059,0.281568627450980,0.287843137254902} error', ...
            '\color[rgb]{0.441568627450980,0.749019607843137,0.432156862745098} correct', ...
            '\color[rgb]{0.600000000000000,0.600000000000000,0.600000000000000} error v correct'};
        
        lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .1;
        set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 6);
    end
end

if doStats
    save(sprintf('%s/Data/GrandAverage/pupilRegressionSignificantCluster.mat', mypath), 'stat');
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


