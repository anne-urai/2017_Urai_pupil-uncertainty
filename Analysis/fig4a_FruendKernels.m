% make all the panels for figure 4a
%     - panel A: frund kernels
%     - panel B: decision strategy for lags 1-3
%     - panel C: valueShift on lag 1-3 with complete stats anova
%     - panel D: frund interaction lags 1-3, explain why we use the multiplicative term
%     - panel D: lag 1 valueShift + Frund interaction bargraph
%     - panel E: decay of neuromodulation over trials (baseline correct the next 7 baselines)
%     - panel F: correct vs error, model comparison

addpath('~/Documents/fieldtrip');
ft_defaults;

clear all; clc; close all;
whichmodulator = 'pupil';
lagGroups = 1:3;

% ========================================================= %
% panel A: Fr?nd kernels for response and stimulus
% ========================================================= %

subplot(341); hold on;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));

a = area(1:3, ones(1, 3) * 0.14, -0.14);
a.FaceColor = [0.98 0.98 0.98];
a.EdgeColor = 'none';
a.LineStyle = 'none';
a.BaseLine.Visible = 'off';

lags = 1:7; subjects = 1:27;
colors = linspecer(9);
colors = colors([5 9], :);
colors(1, :) = colors(1, :) - 0.15; % a bit darker
nlags = 7;
bh = boundedline(lags, nanmean(dat.stimulus), nanstd(dat.stimulus) ./ sqrt(length(subjects)), ...
    lags, nanmean(dat.response), nanstd(dat.response) ./ sqrt(length(subjects)), ...
    'cmap', colors([1 2], :), 'alpha');
xlim([0.5 nlags+0.5]); ylim([-0.15 0.15]); set(gca, 'xtick', lags, 'ytick', -1:0.1:1);
ylabel('History weights'); xlabel('Lags');
text(4, 0.1, 'Response', 'color', colors(2,:), 'fontsize', 7);
text(4.5, -0.1, 'Stimulus', 'color', colors(1,:), 'fontsize', 7);
axis square;

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
cfgstats.numrandomization = 1000; % make sure this is large enough
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
statStim = ft_timelockstatistics(cfgstats, dataStim, dataZero);

%
plot(find(statResp.mask==1), -0.12*ones(1, length(find(statResp.mask==1))), '-', 'color', colors(2, :), 'markersize', 4);
plot(find(statStim.mask==1), -0.13*ones(1, length(find(statStim.mask==1))), '-', 'color', colors(1, :), 'markersize', 4);

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4a_FruendKernels.pdf'));
