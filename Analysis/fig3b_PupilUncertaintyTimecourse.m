function [] = fig3b_PupilUncertaintyTimecourse()
% 3. timecourse of regression betas, separately for correct and error

close all; clear; clc;

% settings
% how to get rid of a possible RT confound?
RTstratification    = true; % 2. add RT as a predictor in designM
doStats             = true; % permutation statistics across the group
plotIndividual      = false; % plots all the individual beta timecourses, takes forever

addpath('~/Documents/fieldtrip');
ft_defaults;

subjects = 1:27;
load('~/Data/pupilUncertainty/GrandAverage/pupilgrandaverage.mat');

warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens
warning('error', 'stats:regress:RankDefDesignMat'); % stop if this happens

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% instead of binning, use regression beta across samples
% 16.08.2015, do this only for the response and feedback parts
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cols = linspecer(4); cols = cols([2 3 1], :); % red green blue

if plotIndividual, figure; end
for sj = unique(subjects),
    
    pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(3).lock.label, 'EyePupil')==1);
    
    % get all timelock
    tl = cat(2, ...
        squeeze(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, :)), ...
        squeeze(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :)));
    
    % do the same for the gazepos, x and y
    xchan    = find(strcmp(pupilgrandavg.timelock{sj}(3).lock.label, 'EyeH')==1);
    gazex = cat(2, ...
        squeeze(pupilgrandavg.timelock{sj}(3).lock.trial(:, xchan, :)), ...
        squeeze(pupilgrandavg.timelock{sj}(4).lock.trial(:, xchan, :)));
    gazex = gazex - nanmean(gazex(:)); % normalize
    
    ychan    = find(strcmp(pupilgrandavg.timelock{sj}(3).lock.label, 'EyeV')==1);
    gazey = cat(2, ...
        squeeze(pupilgrandavg.timelock{sj}(3).lock.trial(:, ychan, :)), ...
        squeeze(pupilgrandavg.timelock{sj}(4).lock.trial(:, ychan, :)));
    gazey = gazey - nanmean(gazey(:)); % normalize
    
    % run a separate glm for correct and error
    cors = [0 1];
    signific = ones(2, size(tl, 2));
    
    for c = 1:2,
        
        thistabledat = pupilgrandavg.timelock{sj}(4).lock.trialinfo;
        trls = find(thistabledat(:, 8) == cors(c));
        
        % add RT as another predictor
        if RTstratification,
            designM = [ones(length(trls), 1) zscore(abs(thistabledat(trls, 4))) ...
                zscore(thistabledat(trls, 6))];
        else % dont include RT
            designM = [ones(length(trls), 1) zscore(abs(thistabledat(trls, 4))) ];
        end

        % regress for each sample
        for s = 1:size(tl, 2),
            
            % add the x and y position of this sample?
            designM2 = [designM gazex(trls, s) gazey(trls, s)];
            
            % when zscoring each sample's timecourse, the nice intercept shape
            % dissappears (duh)
            [b, bint, ~, ~, stats] = regress((tl(trls, s)), designM2);
            grandavg.beta(find(sj==subjects), c, s, :)    = b;
            grandavg.bint(find(sj==subjects), c, s, :)    = [b - bint(:, 1)];
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
save('~/Data/pupilUncertainty/GrandAverage/pupilRegressionBetas.mat', 'grandavg');

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

for whichbeta = 1:size(designM2, 2),
    
    % GRAND average
    sp = subplot(5,3, 1+3*(whichbeta-1));
    % sp = subplot(4,4,3);
    
    % line to indicate zero
    plot([0 size(grandavg.beta, 3)], [0 0], '-', 'color', 'k', 'LineWidth', 0.1);
    
    hold on;
    ph = boundedline(1:size(grandavg.beta, 3), squeeze(nanmean(grandavg.beta(:, :, :, whichbeta)))', ...
        permute(squeeze(nanstd(grandavg.beta(:, :, :, whichbeta))) / sqrt(length(subjects)), [2 3 1]), ...
        'cmap', cols);
    axis tight; ylims = get(gca, 'ylim');
    yval = min(ylims)-0.2*(ylims(2)-ylims(1));
    thisax = gca;
    
    hold on;
    % subfunction to put lines and xlabels at the right spots
    plotLines_respFb(pupilgrandavg.timelock{1}(3).lock, [0 1], ...
        pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
    
    if doStats && whichbeta == 2, % skip stats for the intercept
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
            p(c) = plot(find(stat{c}.mask==1), yval*ones(1, length(find(stat{c}.mask==1))), '.', 'color', cols(c, :), 'markersize', 4);
            yval = yval - 0.08*range(get(thisax, 'ylim'));
            
            % output when this starts to become significant
            signific  = find(stat{c}.mask==1);
            alltiming = [pupilgrandavg.timelock{1}(3).lock.time pupilgrandavg.timelock{1}(4).lock.time];
            disp([alltiming(signific(1)) alltiming(signific(end))]);
        end
        
        % also test their difference
        stat{3} = clusterStat(grandavg_thiscorr(1), grandavg_thiscorr(2), length(subjects));
        p(3) = plot(find(stat{3}.mask==1), yval*ones(1, length(find(stat{3}.mask==1))), '.', 'color', cols(3, :), 'markersize', 4);
        % ylim([-1 0.5]);
    end
    
    % also a shaded area to indicate which part we will use for the
    % statistical comparison
    if 0,
        xticks = get(gca, 'xtick');
        a = area(signific(1):xticks(3), ones(1, length(signific(1):xticks(3))) * max(get(gca, 'ylim')), ...
            min(get(gca, 'ylim')));
        a.FaceColor = [0.8 0.8 0.8];
        a.EdgeColor = 'none';
    end
    set(gca, 'xticklabel', []);
    
    switch whichbeta
        case 1
            set(gca, 'ytick', [-6:2:8]);
            ylabel('\beta intercept');
        case 2
            ylabel('Stimulus difficulty beta');
            set(gca, 'ytick', -1:0.5:0.5);
            % xlabel('Time (ms)');
            
        case 3
            ylabel('\beta reaction time');
            set(gca, 'ytick', -.5:0.5:1);
           % xlabel('Time (ms)');
    end
   % offsetAxes(gca, 0.12, 0);
    
    if whichbeta == 2 && doStats,
        hold on;
        ph2 = plot(180:182, mean(get(thisax, 'ylim'))*ones(3, 10), '.w');
        lh = legend(ph2); % make t
        lh.String = {'\color[rgb]{0.915294117647059,0.281568627450980,0.287843137254902} error', ...
            '\color[rgb]{0.441568627450980,0.749019607843137,0.432156862745098} correct', ...
            '\color[rgb]{0.363921568627451,0.575529411764706,0.748392156862745} error v correct'};
        
        lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .15;
        set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 6);
    end
end

save('~/Data/pupilUncertainty/GrandAverage/pupilRegressionSignificantCluster.mat', 'stat');

% ==================================================================
% show the canonical pupil IRF as well
% ==================================================================

if 0,
    w = 10.1; tmax = 0.930; % from Hoeks and Levelt
    
    % Erlang gamma function from Hoeks and Levelt, Wierda
    t = 0:1/(100):2;
    irf =  (t.^w .* exp(-t.*w ./ tmax));
    
    subplot(6,4,1);
    plot(t, irf); 
    set(gca, 'box', 'off', 'tickdir', 'out', 'yticklabel', []);
    set(gca, 'ytick', get(gca, 'ylim'));
    xlabel('Time (s)'); ylabel({'Impulse'; 'response'});
    offsetAxes;
   % title('Canonical pupil IRF');
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/Fig2c_pupilRegressionTimecourse.pdf'));
disp('SAVED');
end

% ==================================================================
% layout, plots lines to indicate event onset
% ==================================================================
function [] = plotLines(refdata, reftp, stimdata, stimtp, respdata, resptp, fbdata, fbtp)

xticks = []; xlabels = {};
for t = 1:length(reftp),
    xticks = [xticks dsearchn(refdata.time', reftp(t))];
    if reftp(t) == 0,
        xlabels = [xlabels 'Stimulus 1'];
    else
        xlabels = [xlabels reftp(t) * 1000];
    end
end

for t = 1:length(stimtp),
    xticks = [xticks length(refdata.time) + dsearchn(stimdata.time', stimtp(t))];
    if stimtp(t) == 0,
        xlabels = [xlabels 'Stimulus 2'];
    else
        xlabels = [xlabels stimtp(t)* 1000];
    end
end

for t = 1:length(resptp),
    xticks = [xticks length(refdata.time) + length(stimdata.time) + ...
        dsearchn(respdata.time', resptp(t))];
    if resptp(t) == 0,
        xlabels = [xlabels 'Response'];
    else
        xlabels = [xlabels resptp(t)* 1000];
    end
end

for t = 1:length(fbtp),
    xticks = [xticks length(refdata.time) + length(stimdata.time) + ...
        length(respdata.time) + dsearchn(fbdata.time', fbtp(t))];
    if fbtp(t) == 0,
        xlabels = [xlabels 'Feedback'];
    else
        xlabels = [xlabels fbtp(t)* 1000];
    end
end

set(gca, 'XTick', xticks, 'XTickLabel', xlabels, ...
    'XTickLabelRotation', -45, 'tickdir', 'out', 'box', 'off');

% add white lines to indicate transitions between intervals
x = length(refdata.time)+.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);
x = length(refdata.time) + length(stimdata.time) +.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);
x = length(refdata.time) + length(stimdata.time) + length(respdata.time) +.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);

% add dotted  black lines to indicate event onset
x = dsearchn(refdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.1);
x = length(refdata.time) + dsearchn(stimdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.1);
x = length(refdata.time) + length(stimdata.time) + dsearchn(respdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k','LineStyle', '-', 'LineWidth', 0.1);
x = length(refdata.time) + + length(stimdata.time) + ...
    length(respdata.time) + dsearchn(fbdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.1);

end

function [] = plotLines_respFb(respdata, resptp, fbdata, fbtp)

xticks = []; xlabels = {};

for t = 1:length(resptp),
    xticks = [xticks  dsearchn(respdata.time', resptp(t))];
    if resptp(t) == 0,
        xlabels = [xlabels 'Response'];
    else
        xlabels = [xlabels resptp(t)* 1000];
    end
end

for t = 1:length(fbtp),
    xticks = [xticks length(respdata.time) + dsearchn(fbdata.time', fbtp(t))];
    if fbtp(t) == 0,
        xlabels = [xlabels 'Feedback'];
    else
        xlabels = [xlabels fbtp(t)* 1000];
    end
end

set(gca, 'XTick', xticks, 'XTickLabel', xlabels, ...
    'XTickLabelRotation', -45, 'tickdir', 'out', 'box', 'off');

% add white lines to indicate transitions between intervals
x = length(respdata.time) +.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);

% add dotted  black lines to indicate event onset
x = dsearchn(respdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k','LineStyle', '-', 'LineWidth', 0.1);
x =  length(respdata.time) + dsearchn(fbdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.1);

end

function stat = clusterStat(data1, data2, nsubj)

if 1,
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
    cfgstats.numrandomization = 5000; % make sure this is large enough
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
else
    % to do quick plotting tests, dont permute
    for s = 1:size(data1.individual, 3),
        stat.mask(s) = ttest(data1.individual(:, :, s), data2.individual(:, :, s));
    end
end

end


