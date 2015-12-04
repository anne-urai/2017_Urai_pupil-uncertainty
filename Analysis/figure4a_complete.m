% make all the panels for figure 4a
%     - panel A: frund kernels
%     - panel B: decision strategy for lags 1-3
%     - panel C: valueShift on lag 1-3 with complete stats anova
%     - panel D: frund interaction lags 1-3, explain why we use the multiplicative term
%     - panel D: lag 1 valueShift + Frund interaction bargraph
%     - panel E: decay of neuromodulation over trials (baseline correct the next 7 baselines)
%     - panel F: correct vs error, model comparison

clear all; clc; close all;
whichmodulator = 'pupil';
lagGroups = 1:3;

% ========================================================= %
% panel A: Fr?nd kernels for response and stimulus
% ========================================================= %

subplot(341);
load(sprintf('~/Data/UvA_pupil/GrandAverage/historyweights_%s.mat', whichmodulator));

lags = 1:7; subjects = 1:27;
colors = linspecer(8);
nlags = 7;
bh = boundedline(lags, nanmean(dat.stimulus), nanstd(dat.stimulus) ./ sqrt(length(subjects)), ...
    lags, nanmean(dat.response), nanstd(dat.response) ./ sqrt(length(subjects)), ...
    'cmap', colors([2 4], :), 'alpha');
xlim([0.5 nlags+0.5]); ylim([-0.15 0.15]); set(gca, 'xtick', lags, 'ytick', -1:0.1:1);
ylabel('History weights'); xlabel('Lags');
text(4, 0.1, 'Response', 'color', colors(4,:), 'fontsize', 7);
text(4.5, -0.1, 'Stimulus', 'color', colors(2,:), 'fontsize', 7);
plot([1 3], [-.11 -.11], 'k');
axis square;

% ========================================================= %
% panel B: decision strategy for lags 1-3
% ========================================================= %

subplot(342); hold on;

plot([-1 1], [-1 1], 'color', [0.5 0.5 0.5]);
plot([-1 1], [1 -1], 'color', [0.5 0.5 0.5]);

plot(mean(dat.response(:, lagGroups), 2), mean(dat.stimulus(:, lagGroups), 2), ...
    'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
axis tight; axis square; xlims = get(gca, 'xlim'); ylims = get(gca, 'ylim');

maxlim = 0.3;
xlim([-maxlim maxlim]); ylim([-maxlim maxlim]);
xlabel('Response weight'); ylabel('Stimulus weight');
title(sprintf('Lags %d, %d, %d', lagGroups));

% ========================================================= %
% panel C: valueShift on lag 1-3 with complete stats anova
% ========================================================= %

subplot(343);
switch whichmodulator
    case 'pupil'
        whichMod = 'decision_pupil';
    case 'rt'
        whichMod = 'rt';
end

nbins = 3;
spcnt = 0;
subjects = 1:27;
clear grandavg;

for whichLag = 1:7,
    for sj = unique(subjects),
        data = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data2_sj%02d.csv', sj));
        switch whichMod
            case 'rt'
                data.rt = log(data.rt);
        end
        
        % outcome vector need to be 0 1 for logistic regression
        data.resp(data.resp == -1) = 0;
        % remove the trials at the end of each block
        endOfBlockTrls = find(data.trialnr == 50);
        
        % get an overall logistic fit
        [b, dev, stats] = glmfit(data.motionstrength, data.resp, ...
            'binomial','link','logit');
        grandavg.overallLogistic(sj, :) = b;
        
        % previous response
        resps = [0 1];
        for r = 1:2,
            
            if nbins > 2,
                uncQs = quantile(data.(whichMod)(data.resp == resps(r)), nbins - 1);
            elseif nbins == 2,
                uncQs = median(data.(whichMod)(data.resp == resps(r)));
            end
            
            % uncertainty bins
            for u = 1:nbins,
                
                switch u
                    case 1
                        trls = find(data.resp == resps(r) & data.(whichMod) < uncQs(u));
                    case nbins
                        trls = find(data.resp == resps(r) & data.(whichMod) > uncQs(u-1));
                    otherwise
                        trls = find(data.resp == resps(r) & ...
                            data.(whichMod) > uncQs(u-1) & data.(whichMod) < uncQs(u));
                end
                
                % only use correct trials!
                % trls = intersect(trls, find(data.correct == 1));
                
                % with this selection, take the trials after that
                laggedtrls = trls+whichLag;
                % exclude trials at the end of the block
                if any(laggedtrls > size(data, 1)),
                    trls(laggedtrls > size(data, 1)) = [];
                    laggedtrls(laggedtrls > size(data, 1)) = [];
                end
                
                % remove trials that dont match in block nr
                removeTrls = data.blocknr(laggedtrls) ~= data.blocknr(trls);
                laggedtrls(removeTrls) = [];
                trls(removeTrls) = [];
                
                % fit logistic regression
                thisdat = data(laggedtrls, :);
                
                [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                    'binomial','link','logit');
                
                % save betas
                grandavg.logistic(sj, r, u, whichLag, :) = b;
                
            end % uncertainty bin
        end % resp
    end % sj
end

% normalize by their overall bias
grandavg.logistic(:, :, :, :, 1) = bsxfun(@minus, grandavg.logistic(:, :, :, :, 1), grandavg.overallLogistic(:, 1));
% take the mean over lags
grandavg.groupedLogistic = squeeze(nanmean(grandavg.logistic(:, :, :, lagGroups, :), 4)); % average over lags

colors = linspecer(6);
colors = colors([5 6], :);
stimx2 = 1:u;

subplot(3,4,3);
hold on;
plot([1 max(stimx2)], [0 0], '-', 'color', [0.4 0.4 0.4], 'LineWidth', 0.1);

for r = [1 2],
    h = ploterr(stimx2, ...
        squeeze(nanmean(grandavg.groupedLogistic(:, r, :, 1))),  [], ...
        squeeze(nanstd(grandavg.groupedLogistic(:, r, :, 1))) ./ sqrt(length(subjects)), ...
        '-',  'hhxy', 0.001);
    set(h(1), 'color', colors(r, :), ...
        'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
    set(h(2), 'color', colors(r, :));
    stimx2 = stimx2 - 0.15;
end

ylabel({'Response bias_{t0}'});
set(gca, 'ytick', [-.15 0  0.15]);
whichModTitle = whichMod; whichModTitle(1) = upper(whichModTitle(1));
xlim([0.5 3.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});
ylim([-.15 .15]);
switch whichmodulator
    case 'rt'
        xlabel('Reaction time_{t-1:-3}');
    otherwise
        xlabel('Decision pupil_{t-1:-3}');
end

% repeated measures anova on those bins
cnt = 0; clear x s f
for sj = 1:27,
    for r = 1:2,
        for u = 1:3,
            cnt = cnt + 1;
            x(cnt) = grandavg.groupedLogistic(sj, r, u, 1);
            s(cnt) = sj;
            f{1}(cnt) = r; % factor of previous response
            f{2}(cnt) = u; % factor of previous uncertainty
        end
    end
end

stats = rm_anova(x, s, f);

subplot(3,4,4);
text(0, 0.9, sprintf('Main RespBias F_{(%d,%d)} = %.2f', stats.f1.df(1), stats.f1.df(2), stats.f1.fstats));
text(0, 0.75,  sprintf('p = %.3f', stats.f1.pvalue));

text(0, 0.5, sprintf('Main modul F_{(%d,%d)} = %.2f', stats.f2.df(1), stats.f2.df(2), stats.f2.fstats));
text(0, 0.35,  sprintf('p = %.3f', stats.f2.pvalue));

text(0, 0.2, sprintf('Interaction F_{(%d,%d)} = %.2f', stats.f1xf2.df(1), stats.f1xf2.df(2), stats.f1xf2.fstats));
text(0, 0.05,  sprintf('p = %.3f', stats.f1xf2.pvalue));
axis off;

% ========================================================= %
% panel D: Fruend weights for resp and pupil lag 1-3 grouped
% ========================================================= %

subplot(3,4,5); hold on;

% response weight
bar(1, squeeze(mean(mean(dat.response(:, lagGroups), 2))), ...
    'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none', 'barwidth', 0.5);
h = ploterr(1, squeeze(mean(mean(dat.response(:, lagGroups), 2))), ...
    [], squeeze(std(mean(dat.response(:, lagGroups), 2))) / sqrt(length(subjects)), ...
    'k.', 'hhxy', 0.0000000000000001);
set(h(1), 'marker', 'none');
axis tight; axis square;

% stats!
[pval] = randtest1d(mean(dat.response(:, lagGroups), 2));
sigstar({[1 1]}, pval);

% start with a graph for all the trials combined
bar(2, squeeze(mean(mean(dat.response_pupil(:, lagGroups), 2))), ...
    'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none', 'barwidth', 0.5);
h = ploterr(2, squeeze(mean(mean(dat.response_pupil(:, lagGroups), 2))), ...
    [], squeeze(std(mean(dat.response_pupil(:, lagGroups), 2))) / sqrt(length(subjects)), ...
    'k.', 'hhxy', 0.0000000000000001);
set(h(1), 'marker', 'none');
axis tight; axis square;

[pval] = randtest1d(mean(dat.response_pupil(:, lagGroups), 2));
sigstar({[2 2]}, pval);

switch whichmodulator
    case 'rt'
        set(gca, 'xtick', 1:2, 'xticklabel', {'resp', 'resp*rt'});
    case 'pupil'
        set(gca, 'xtick', 1:2, 'xticklabel', {'resp', 'resp*pupil'});
end
xlabel('Lags 1-3');
ylabel('Betas Fruend');
ylim([-0.05 0.1]); xlim([0.5 2.5]);

% ========================================================= %
% panel E: neuromodulatory decay
% ========================================================= %

subjects = 1:27;
clear neuromodDecay
for sj = unique(subjects),
    data = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data2_sj%02d.csv', sj));
    trls = 1:size(data, 1);
    
    for lag = 1:10,
        
        % with this selection, take the trials after that
        laggedtrls = trls+lag;
        
        % exclude trials at the end of the block
        if any(laggedtrls > size(data, 1)),
            trls(laggedtrls > size(data, 1)) = [];
            laggedtrls(laggedtrls > size(data, 1)) = [];
        end
        % remove trials that dont match in block nr
        noblockmatch = data.blocknr(laggedtrls) ~= data.blocknr(trls);
        laggedtrls(noblockmatch) = [];
        trls(noblockmatch) = [];
        
        % subtract correct the baseline on tLagged by the baseline on t0
        baselineCorr = projectout(data.baseline_pupil(laggedtrls), data.baseline_pupil(trls));
        % subtract or project out?
        baselineCorr = data.baseline_pupil(laggedtrls) - data.baseline_pupil(trls);
        
        % these look more or less the same
        neuromodDecay.spear(sj, lag) = corr(baselineCorr, data.decision_pupil(trls), 'type', 'Pearson');
        neuromodDecay.pear(sj, lag) = corr(baselineCorr, data.decision_pupil(trls), 'type', 'Spearman');
        
    end
end

subplot(3,4,6); hold on;
b1 = boundedline(1:lag, mean(neuromodDecay.pear), std(neuromodDecay.pear) ./ sqrt(length(subjects)), 'cmap', [0.8 0.8 0.8]);
b2 = boundedline(1:lag, mean(neuromodDecay.spear), std(neuromodDecay.spear) ./ sqrt(length(subjects)), 'cmap', [0.3 0.3 0.3]);
axis tight;  xlim([0.5 lag+0.5]);
xlabel('Lags'); ylabel({'Correlation'});

lh = legend([b1 b2], {'Pearson', 'Spearman'}, 'Location', 'NorthEast'); legend boxoff;

% ========================================================= %
% panel F: pupil valueShift for lag 1
% ========================================================= %

% take the mean over lags
lagGroups = 1;
grandavg.groupedLogistic = squeeze(nanmean(grandavg.logistic(:, :, :, lagGroups, :), 4)); % average over lags

colors = linspecer(6);
colors = colors([5 6], :);
stimx2 = 1:u;

subplot(347);
hold on;
plot([1 max(stimx2)], [0 0], '-', 'color', [0.4 0.4 0.4], 'LineWidth', 0.1);

for r = [1 2],
    h = ploterr(stimx2, ...
        squeeze(nanmean(grandavg.groupedLogistic(:, r, :, 1))),  [], ...
        squeeze(nanstd(grandavg.groupedLogistic(:, r, :, 1))) ./ sqrt(length(subjects)), ...
        '-',  'hhxy', 0.001);
    set(h(1), 'color', colors(r, :), ...
        'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
    set(h(2), 'color', colors(r, :));
    stimx2 = stimx2 - 0.15;
end

ylabel({'Response bias_{t0}'});
set(gca, 'ytick', [-.15 0  0.15]);
whichModTitle = whichMod; whichModTitle(1) = upper(whichModTitle(1));
xlim([0.5 3.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});
ylim([-.15 .15]);
switch whichmodulator
    case 'rt'
        xlabel('Reaction time_{t-1}');
    otherwise
        xlabel('Decision pupil_{t-1}');
end

%% repeated measures anova on those bins
cnt = 0; clear x s f
for sj = 1:27,
    for r = 1:2,
        for u = 1:3,
            cnt = cnt + 1;
            x(cnt) = grandavg.groupedLogistic(sj, r, u, 1);
            s(cnt) = sj;
            f{1}(cnt) = r; % factor of previous response
            f{2}(cnt) = u; % factor of previous uncertainty
        end
    end
end

stats = rm_anova(x, s, f);

subplot(3,4,8);
text(0, 0.9, sprintf('Main RespBias F_{(%d,%d)} = %.2f', stats.f1.df(1), stats.f1.df(2), stats.f1.fstats));
text(0, 0.75,  sprintf('p = %.3f', stats.f1.pvalue));

text(0, 0.5, sprintf('Main modul F_{(%d,%d)} = %.2f', stats.f2.df(1), stats.f2.df(2), stats.f2.fstats));
text(0, 0.35,  sprintf('p = %.3f', stats.f2.pvalue));

text(0, 0.2, sprintf('Interaction F_{(%d,%d)} = %.2f', stats.f1xf2.df(1), stats.f1xf2.df(2), stats.f1xf2.fstats));
text(0, 0.05,  sprintf('p = %.3f', stats.f1xf2.pvalue));
axis off;

% ========================================================= %
% panel H: lag 1 model values
% ========================================================= %

clc; clearvars -except whichmodulator;

subjects = 1:27;
grandavg.betas = nan(27, 3, 7);
grandavg.pval = nan(27, 3, 7);

for sj = unique(subjects),
    data = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data2_sj%02d.csv', sj));
    
    % get the variables we need
    resp = data.resp; resp(resp == -1) = 0; % predict response identity
    motionstrength = zscore(data.motionstrength);
    
    switch whichmodulator
        case 'pupil'
            prevPupil      = circshift(zscore(data.decision_pupil), 1);
            xlabs = {'bias', 'stim', 'resp','resp*pupil'};
        case 'rt'
            prevPupil      = circshift(zscore(log(data.rt)), 1);
            xlabs = {'bias', 'stim', 'resp','resp*rt'};
    end
    prevResp       = circshift(data.resp, 1);
    prevStim       = circshift(data.stim, 1);
    
    % make design matrix - intercept will be added automatically
    designM = [motionstrength prevResp prevResp.*prevPupil];
    
    % don't use trials that are at the beginning of each block
    trlDif = [0; diff(data.trialnr)];
    
    removeTrls = false(size(trlDif));
    removeTrls(trlDif < 1) = true;
    removeTrls(trlDif > 1) = true;
    removeTrls(find(trlDif > 1) + 1) = true;
    
    % in these trials, the history wont be able to predict the response
    designM(removeTrls==1, 2:end) = 0;
    
    % do this separately for error, correct and all trials together
    cors = [0 1];
    for c = 1:3,
        correctPrev = circshift(data.correct, 1);
        if c == 3,
            trls = 1:length(resp);
        else
            trls = find(correctPrev == cors(c));
        end
        
        mdl  = fitglm(designM(trls, :), resp(trls), ...
            'distr', 'binomial', 'link', 'logit');
        grandavg.betas(sj, c, 1:length(mdl.Coefficients.Estimate)) = mdl.Coefficients.Estimate;
        grandavg.pval(sj, c, 1:length(mdl.Coefficients.pValue))  = mdl.Coefficients.pValue;
        
        grandavg.dev(sj, 3) = mdl.Deviance;
        
        if c == 3, % only on all trials
            
            % get the deviance for the null model without history
            mdl  = fitglm(designM(trls, 1), resp(trls), ...
                'distr', 'binomial', 'link', 'logit');
            grandavg.dev(sj, 1) = mdl.Deviance;
            
            % model with  history
            mdl  = fitglm(designM(trls, 1:2), resp(trls), ...
                'distr', 'binomial', 'link', 'logit');
            grandavg.dev(sj, 2) = mdl.Deviance;
            
        end
    end
end

subplot(3,4,9); hold on;
colors = linspecer(3); colors = colors([2 3 1], :);
what2plot = [3 4];

% start with a graph for all the trials combined
bar(1:length(what2plot), squeeze(mean(grandavg.betas(:, 3, what2plot))), ...
    'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none', 'barwidth', 0.5);
h = ploterr(1:length(what2plot), squeeze(mean(grandavg.betas(:, 3, what2plot))), ...
    [], squeeze(std(grandavg.betas(:, 3, what2plot))) / sqrt(length(subjects)), ...
    'k.', 'hhxy', 0.0000000000000001);
set(h(1), 'marker', 'none');
axis tight; axis square;
set(gca, 'xlim', [0.5 2.5]);
set(gca, 'xcolor', 'k', 'ycolor', 'k');
set(gca, 'xtick', [1 2], 'xticklabel', xlabs(what2plot));
xlabel('Lag 1');

clear pval;
for n = 1:length(what2plot),
    %     pval(n) = randtest1d(, ...
    %         zeros(size(squeeze(grandavg.betas(:, c, what2plot(n))))), ...
    %         0,10000);
    [pval(n)] = randtest1d(squeeze(grandavg.betas(:, c, what2plot(n))));
end

% Odds ratios and beta coefficients both estimate the effect of an exposure on the
% outcome, the later one being the natural logarithm of the former one.
% beta = log(OR); - should be normally distributed!

sigstar({[1 1], [2 2]}, pval);
ylabel('Betas LogRes');
box off;

% ========================================================= %
% panel I: lag 1 model values correct vs error
% ========================================================= %

% now, do the same but separately for correct and error
subplot(3,4,10); hold on;
clear pval;
for sp = 2,
    
    if sp == 1,
        what2plot = 1:2; % bias stimswitch
    elseif sp == 2,
        what2plot = 3:4; % resp pupil interaction
    end
    
    % correct and error
    for c = 1:2,
        
        switch c
            case 1
                x = 1:2:length(what2plot)*2;
            case 2;
                x = 2:2:length(what2plot)*2;
        end
        
        bar(x, squeeze(mean(grandavg.betas(:, c, what2plot))), ...
            'facecolor', colors(c, :), 'edgecolor', 'none', 'barwidth', 0.4);
        h = ploterr(x, squeeze(mean(grandavg.betas(:, c, what2plot))), ...
            [], squeeze(std(grandavg.betas(:, c, what2plot))) / sqrt(length(subjects)), ...
            'k.', 'hhxy', 0.0000000000000001);
        set(h(1), 'marker', 'none');
        set(gca, 'xlim', [0.5 4.5]);
        
        % within error and correct, compute significance across the group
        for n = 1:length(what2plot),
            [pval(x(n))] = randtest1d(squeeze(grandavg.betas(:, c, what2plot(n))));
        end
    end
    
    % also add sigstars to compare each regressor between error and correct
    pidx = [5 6];
    for n = 1:length(what2plot),
        [pval(pidx(n))] = randtest1d(squeeze(grandavg.betas(:, 1, what2plot(n))), squeeze(grandavg.betas(:, 2, what2plot(n))));
    end
    sigstar({[1 1], [2 2], [3 3], [4 4], [1 2], [3 4]}, pval, 1);
    
    set(gca, 'xlim', [0.5 4.5]);
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    set(gca, 'xtick', 1.5:2:length(what2plot)*2, 'xticklabel', xlabs(what2plot));
    xlabel('Lag 1');
    
    ylabel('Betas LogRes');
    box off;
end

% ========================================================= %
% between subject scatter
% ========================================================= %
load(sprintf('~/Data/UvA_pupil/GrandAverage/historyweights_%s.mat', whichmodulator));

mainW = grandavg.betas(:, 3, 3);
intW  = grandavg.betas(:, 3, 4);
colors = linspecer(3);

subplot(3,4,11);
plot(mainW, intW, ...
    'o', 'MarkerFaceColor', colors(3,:) , 'MarkerEdgeColor', 'w');
box off; xlabel('Correct response');
switch whichmodulator
    case 'rt';
        ylabel('Correct resp*rt');
    case 'pupil'
        ylabel('Correct resp*pupil');
end
[rho, pval] = corr(mainW, intW, 'type', 'p');
if pval < 0.05,
    ls = lsline; set(ls, 'color', [0.7 0.7 0.7]);
    %  assert(1==0)
    disp([rho pval]);
    text(0.3, 0, sprintf('rho = %.3f', rho));
    text(0.3, -0.05, sprintf('p = %.3f', pval));
end
title('LogRes weights');
axis square;

mainW = dat.correct(:, 1);
intW  = dat.correct_pupil(:, 1);
colors = linspecer(3);

subplot(3,4,12);
plot(mainW, intW, ...
    'o', 'MarkerFaceColor', colors(3,:) , 'MarkerEdgeColor', 'w');
box off; xlabel('Correct response');
switch whichmodulator
    case 'rt';
        ylabel('Correct resp*rt');
    case 'pupil'
        ylabel('Correct resp*pupil');
end
[rho, pval] = corr(mainW, intW, 'type', 'pearson');
if pval < 0.05,
    ls = lsline; set(ls, 'color', [0.7 0.7 0.7]);
    %  assert(1==0)
    disp([rho pval]);
    text(0.3, 0, sprintf('rho = %.3f', rho));
    text(0.3, -0.05, sprintf('p = %.3f', pval));
end
title('Fruend weights');
axis square;

%% ========================================================= %
% save
% ========================================================= %

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/completefigure4_%s.pdf', whichmodulator));
