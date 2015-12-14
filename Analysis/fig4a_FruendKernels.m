function [] = fig4a_FruendKernels(lagGroups, whichmodulator, grouping)

% make all the panels for figure 4a
%     - panel A: frund kernels
%     - panel B: decision strategy for lags 1-3
%     - panel C: valueShift on lag 1-3 with complete stats anova
%     - panel D: frund interaction lags 1-3, explain why we use the multiplicative term
%     - panel D: lag 1 valueShift + Frund interaction bargraph
%     - panel E: decay of neuromodulation over trials (baseline correct the next 7 baselines)
%     - panel F: correct vs error, model comparison

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end

addpath('~/Documents/fieldtrip');
ft_defaults;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

switch grouping
    case 'individual'
        %subplot(441);
        hold on;
        colors = linspecer(9);
        posRespSj = find(dat.response(:, 1) > 0);
        negRespSj = find(dat.response(:, 1) < 0);
        
        plot(dat.response(posRespSj, :)', 'color', colors(8, :));
        plot(ones(size(dat.response(posRespSj))), dat.response(posRespSj, 1), ...
            'o', 'markerfacecolor', colors(8,:), 'markeredgecolor', 'w', 'markersize', 4);
        plot(dat.response(negRespSj, :)', 'color', colors(9, :));
        plot(ones(size(dat.response(negRespSj))), dat.response(negRespSj, 1), ...
            'o', 'markerfacecolor', colors(9,:), 'markeredgecolor', 'w', 'markersize', 4);
        xlim([0.5 7]); ylim([-.4 .4]); set(gca, 'xtick', 1:7, 'ytick', -1:0.2:1);
        ylabel('Response weight'); xlabel('Lags');
        axis square;
        
        % ========================================================= %
        % panel A: Fr?nd kernels for response and stimulus
        % ========================================================= %
        
    case 'avg'
        %subplot(442); 
        hold on;
        
        if lagGroups == 1,
            lagGroups = [0.8 1.2];
        end
        
        a = area(lagGroups, ones(1, length(lagGroups)) * 0.14, -0.01);
        a.FaceColor = [0.9 0.9 0.9];
        a.EdgeColor = 'none';
        a.LineStyle = 'none';
        a.BaseLine.Visible = 'off';
        
        lags = 1:7; subjects = 1:27;
        colors = linspecer(9);
        colors = colors([2], :);
        
        colors(1, :) = colors(1, :) - 0.15; % a bit darker
        nlags = 7;
        bh = boundedline(lags, nanmean(dat.response), nanstd(dat.response) ./ sqrt(length(subjects)), ...
            'cmap', colors, 'alpha');
        xlim([0.5 nlags]); ylim([-0.02 0.15]); set(gca, 'xtick', lags, 'ytick', -1:0.1:1);
        ylabel('Response weights'); %xlabel('Lags');
        %text(4, 0.1, 'Response', 'color', colors(2,:), 'fontsize', 7);
        %text(4.5, -0.1, 'Stimulus', 'color', colors(1,:), 'fontsize', 7);
        axis square;
        
        
        if 0,
            % cluster based permutation test and indicate stats.mask
            cfgstats                  = [];
            cfgstats.method           = 'montecarlo'; % permutation test
            cfgstats.statistic        = 'ft_statfun_depsamplesT'; % dependent samples ttest
            
            % do cluster correction
            cfgstats.correctm         = 'cluster';
            cfgstats.clusteralpha     = 0.05;
            cfgstats.tail             = 0; % two-tailed!
            cfgstats.clustertail      = 0; % two-tailed!
            cfgstats.alpha            = 0.025;
            cfgstats.numrandomization = 100; % make sure this is large enough
            cfgstats.randomseed       = 1; % make the stats reproducible!
            
            % use only our preselected sensors for the time being
            cfgstats.channel          = 'weights';
            
            % specifies with which sensors other sensors can form clusters
            cfgstats.neighbours       = []; % only cluster over data and time
            
            design = zeros(2,2*length(subjects));
            for i = 1:length(subjects),  design(1,i) = i;        end
            for i = 1:length(subjects),  design(1,length(subjects)+i) = i;  end
            design(2,1:length(subjects))         = 1;
            design(2,length(subjects)+1:2*length(subjects)) = 2;
            
            cfgstats.design   = design;
            cfgstats.uvar     = 1;
            cfgstats.ivar     = 2;
            
            % make resp data
            dataResp.time       = 1:7;
            dataResp.label      = {'weights'};
            dataResp.dimord     = 'subj_chan_time';
            dataResp.individual = permute(dat.response, [1 3 2]);
            
            % make stim data
            dataStim = dataResp;
            dataStim.individual = permute(dat.stimulus, [1 3 2]);
            
            % compare against null
            dataZero = dataResp;
            dataZero.individual = zeros(size(dataResp.individual));
            
            statResp = ft_timelockstatistics(cfgstats, dataResp, dataZero);
            %statStim = ft_timelockstatistics(cfgstats, dataStim, dataZero);
            
            %
            plot(find(statResp.mask==1), -0.01*ones(1, length(find(statResp.mask==1))), '-', 'color', colors(1, :), 'markersize', 4);
            %plot(find(statStim.mask==1), -0.13*ones(1, length(find(statStim.mask==1))), '-', 'color', colors(1, :), 'markersize', 4);
            
        else % dont do this every time for speed
            plot(2:7, -0.01*ones(1, 6), '-', 'color', colors(1, :), 'markersize', 4);
        end
end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4a_FruendKernels.pdf'));
