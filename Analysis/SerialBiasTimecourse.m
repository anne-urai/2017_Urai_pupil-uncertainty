function [] = SerialBiasTimecourse(plotAll)
% timecourse of predictive information

global mypath;
if ~exist('plotall', 'var'); plotAll = false; end % in principle, only the main evidence strength regressor

subjects = 1:27;
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens
warning('error', 'stats:regress:RankDefDesignMat'); % stop if this happens

% ========================================================
% instead of binning, use regression beta across samples
% ========================================================

colors = cbrewer('qual', 'Set1', 9);
colors = colors([1 2 9], :); % red blue grey

if exist(sprintf('%s/Data/GrandAverage/SerialBiasRegressionBetas.mat', mypath), 'file'), ...
        load(sprintf('%s/Data/GrandAverage/SerialBiasRegressionBetas.mat', mypath)); % dont redo, takes time;
else
    for sj = unique(subjects),
        disp(sj);
        pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(3).lock.label, 'EyePupil')==1);
        
        % get all timelock
        tl = cat(2, ...
            squeeze(pupilgrandavg.timelock{sj}(1).lock.trial(:, pupilchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(2).lock.trial(:, pupilchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, :)), ...
            squeeze(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :)));
        
        % make design matrix
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        
        % get the variables we need
        data.motionstrength = zscore(data.motionstrength);
        data.prevResp       = circshift(data.resp, 1);
        data.resp(data.resp == -1) = 0; % predict response identity
        data.prevRT         = circshift(zscore(data.rtNorm), 1);
        
        % don't use trials that are at the beginning of each block
        trlDif = [0; diff(data.trialnr)];
        removeTrls = false(size(trlDif));
        removeTrls(trlDif < 1) = true;
        removeTrls(trlDif > 1) = true;
        removeTrls(find(trlDif > 1) + 1) = true;
        data.resp(removeTrls)        = NaN;
        
        nanzscore = @(x) (x - nanmean(x)) ./ nanstd(x);
        
        % regress for each sample
        for s = 1:size(tl, 2),
            
            % get previous pupil at this one sample
            data.prevPupil      = circshift(nanzscore(tl(:, s)), 1);
            
            % run regression
            mdl = fitglm(data, 'resp ~ motionstrength + prevResp*prevPupil + prevResp*prevRT', ...
                'distr', 'binomial', 'link', 'logit');
            
            grandavg.betas(sj, s, :) = mdl.Coefficients.Estimate;
        end
    end
    
    % save for next time
    grandavg.names = mdl.CoefficientNames;
    savefast(sprintf('%s/Data/GrandAverage/SerialBiasRegressionBetas.mat', mypath), 'grandavg'); % dont redo, takes time;
end

% ========================================================
% plot the timecourse of regression coefficients
% ========================================================

whichBetas2plot = [7]; % prevResp*prevPupil and prevResp*prevPupil
colors = [0 0 0];

for b = 1:length(whichBetas2plot),
    thisb = whichBetas2plot(b);
    subplot(2,2,[1:2] + (b-1)*2);
    
    % line to indicate zero
    plot([0 size(grandavg.betas, 2)], [0 0], ':', 'color', 'k', 'LineWidth', 0.2);
    
    hold on;
    boundedline(1:size(grandavg.betas, 2), squeeze(nanmean(grandavg.betas(:, :, thisb)))', ...
        permute(squeeze(nanstd(grandavg.betas(:, :, thisb))) / sqrt(length(subjects)), [2 3 1]), ...
        'cmap', colors);
    axisNotSoTight;
    
    hold on;
    % subfunction to put lines and xlabels at the right spots
    plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
        pupilgrandavg.timelock{1}(2).lock, [0], ...
        pupilgrandavg.timelock{1}(3).lock, [0 1], ...
        pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
    
    
    if 0,
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
        
        % plot
        % plot on top
        p(c) = plot(find(stat{c}.mask==1), yval*ones(1, length(find(stat{c}.mask==1))), '.', 'color', colors(c, :), 'markersize', 4);
        yval = yval - 0.08*range(get(thisax, 'ylim'));
        
        % output when this starts to become significant
        signific  = find(stat{c}.mask==1);
        alltiming = [pupilgrandavg.timelock{1}(3).lock.time pupilgrandavg.timelock{1}(4).lock.time];
        try
            disp([alltiming(signific(1)) alltiming(signific(end))]);
        end
        p(3) = plot(find(stat{3}.mask==1), yval*ones(1, length(find(stat{3}.mask==1))), '.', 'color', colors(3, :), 'markersize', 4);
    end
    
    ylabel('Beta weight (a.u.)');
    title(grandavg.names{thisb});
    print(gcf, '-dpdf', sprintf('%s/Figures/SerialBiasTimecourse.pdf', mypath));

    %
    %     hold on;
    %     ph2 = plot(180:182, mean(get(thisax, 'ylim'))*ones(3, 10), '.w');
    %     lh = legend(ph2);
    %     legnames = {'error', 'correct', 'error v correct'};
    %     for i = 1:length(legnames),
    %         str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
    %     end
    %     lh.String = str;
    %     lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .1;
    %     lpos(2) = lpos(2) - 0.08;
    %     set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 6);
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


