%% 1. dependence on pupil-based uncertainty, all trials
% does this pattern increase with higher uncertainty on previous trial?
clear; clc; close;
nbins = 3;
subjects = 1:27;
whichLag = 1;
whichMod = 'baseline_pupil';
spcnt = 0;

spcnt = spcnt + 1;
for sj = unique(subjects),
    
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    switch whichMod
        case 'rt'
            data.rt = log(data.rt);
        case 'baseline_pupil'
            data.baseline_pupil = circshift(data.baseline_pupil, -1);
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
            usetrls = trls+whichLag;
            % exclude trials at the end of the block
            if any(usetrls > size(data, 1)),
                trls(usetrls > size(data, 1)) = [];
                usetrls(usetrls > size(data, 1)) = [];
            end
            
            % remove trials that dont match in block nr
            usetrls(data.blocknr(usetrls) ~= data.blocknr(trls)) = [];
            
            % fit logistic regression
            thisdat = data(usetrls, :);
            
            [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                'binomial','link','logit');
            stimx = -4:0.01:4;
            yfit = glmval([b(1) - grandavg.overallLogistic(sj, 1); b(2)], stimx, 'logit');
            %  plot(x, yfit);
            
            % save betas
            grandavg.logistic(sj, r, u, :) = b;
            grandavg.curve(sj, r, u, :)    = yfit; % curve
            
        end % uncertainty bin
    end % resp
end % sj

%% normalize by their overall bias
grandavg.logistic(:, :, :, 1) = bsxfun(@minus, grandavg.logistic(:, :, :, 1), grandavg.overallLogistic(:, 1));

%% plot 2. plot interaction effect and 2way anova
colors = linspecer(4);
colors = colors([1 4], :);
stimx2 = 1:u;

subplot(4,4,spcnt); spcnt = spcnt + 1;
hold on;
plot([1 max(stimx2)], [0 0], '-', 'color', [0.4 0.4 0.4], 'LineWidth', 0.1);

for r = [1 2],
    h = ploterr(stimx2, ...
        squeeze(nanmean(grandavg.logistic(:, r, :, 1))),  [], ...
        squeeze(nanstd(grandavg.logistic(:, r, :, 1))) ./ sqrt(length(subjects)), ...
        '-',  'hhxy', 0.001);
    set(h(1), 'color', colors(r, :), ...
        'marker', '.', 'markerfacecolor', colors(r, :), 'markersize', 12);
    set(h(2), 'color', colors(r, :));
    stimx2 = stimx2 - 0.15;
end

ylabel({'Response bias_{t0}'});
set(gca, 'ytick', [-.15 0  0.15]);
whichModTitle = whichMod; whichModTitle(1) = upper(whichModTitle(1));
xlabel({'Baseline pupil_{0}'});
axis square;
xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});
text(2, -.12, sprintf('Choice_{t-%d} B', whichLag), 'color', colors(1, :));
text(2, .12, sprintf('Choice_{t-%d} A', whichLag), 'color', colors(2, :));
ylim([-.2 .2]);

%% repeated measures anova on those bins
cnt = 0; clear x s f
for sj = 1:27,
    for r = 1:2,
        for u = 1:nbins,
            cnt = cnt + 1;
            x(cnt) = grandavg.logistic(sj, r, u, 1);
            s(cnt) = sj;
            f{1}(cnt) = r; % factor of previous response
            f{2}(cnt) = u; % factor of previous uncertainty
        end
    end
end

stats = rm_anova(x, s, f);

%    >> stats = rm_anova(x, s, f);
%           x = 1xN vector with data values
%           s = 1xN vector with subject numbers, or the factor to repeat
%           measures over
%           f = factors, 1xNrFactors cell array. Each cell array should
%           have an 1xN vector with the level within this factor.

disp(stats.f1xf2)
% title({sprintf('F_{(%d,%d)} = %.2f', stats.f1xf2.df(1), stats.f1xf2.df(2), stats.f1xf2.fstats); ...
%    sprintf('p = %.3f', stats.f1xf2.pvalue)});

%% plot 3. plot main effect on response repetition and one-way anova

resp1 = -squeeze(grandavg.logistic(:, 1, :, 1));
resp2 = squeeze(grandavg.logistic(:, 2, :, 1));
grandavg.logisticRep = (resp1 + resp2) ./ 2;

colors = linspecer(5);
stimx2 = 1:u;

subplot(4,4,spcnt); hold on; spcnt = spcnt + 1;
h = ploterr(stimx2, ...
    squeeze(nanmean(grandavg.logisticRep)),  [], ...
    squeeze(nanstd(grandavg.logisticRep)) ./ sqrt(length(subjects)), ...
    '-',  'hhxy', 0.001);
set(h(1), 'color', colors(5, :), ...
    'marker', '.', 'markerfacecolor', colors(5, :), 'markersize', 12);
set(h(2), 'color', colors(5, :));
plot([1 max(stimx2)], [0 0], '-', 'color', [0.4 0.4 0.4], 'LineWidth', 0.1);

ylabel({'Repetition bias_{0}'});
ylim([-.07 .15]);
xlabel({'Baseline pupil_{0}'});
axis square;
xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:3, 'xticklabel', {'low', 'med', 'high'});

if 0,
    %% repeated measures anova on those bins
    cnt = 0; clear x s f
    for sj = 1:27,
        for u = 1:3,
            cnt = cnt + 1;
            x(cnt) = grandavg.logisticRep(sj, u);
            s(cnt) = sj;
            f{1}(cnt) = u; % factor of previous uncertainty
        end
    end
    clear stats
    stats = rm_anova(x, s, f);
    
    %    >> stats = rm_anova(x, s, f);
    %           x = 1xN vector with data values
    %           s = 1xN vector with subject numbers, or the factor to repeat
    %           measures over
    %           f = factors, 1xNrFactors cell array. Each cell array should
    %           have an 1xN vector with the level within this factor.
    
end

% disp(stats.f1)
subplot(4,4,spcnt);
text(0.2, 0.8, sprintf('Main RespB F_{(%d,%d)} = %.2f', stats.f1.df(1), stats.f1.df(2), stats.f1.fstats));
text(0.2, 0.6,  sprintf('p = %.3f', stats.f1.pvalue));

text(0.2, 0.3, sprintf('Interaction F_{(%d,%d)} = %.2f', stats.f1xf2.df(1), stats.f1xf2.df(2), stats.f1xf2.fstats));
text(0.2, 0.1,  sprintf('p = %.3f', stats.f1xf2.pvalue));
axis off;


%% save fig
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/FigS_psychFuncShift_pupilBaseline_%s.pdf', whichMod));
