function fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, grouping)
% plot both the effect of pupil on overall repetition bias and show that
% this is symmetrical for both previous choices

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('grouping', 'var'); grouping = 'all'; end

switch whichmodulator
    case 'pupil'
        whichMod = 'decision_pupil';
    case 'feedbackpupil'
        whichMod = 'feedback_pupil';
    case 'rt'
        whichMod = 'rt';
    case 'fb-decpupil'
        % see below for projecting out
        whichMod = 'feedback_pupil';
end

nbins = 3;
subjects = 1:27;
whichLags = 1:3; % lag 1

for lag = whichLags,
    for sj = unique(subjects),
        data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
        
        % in this case, take out decision effects
        switch whichmodulator
            case 'fb-decpupil'
                data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        end
        
        % outcome vector need to be 0 1 for logistic regression
        data.resp(data.resp == -1) = 0;
        
        % get an overall logistic fit
        grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, data.resp, ...
            'binomial','link','logit');
        
        % previous response
        if nbins > 2,
            uncQs = quantile(data.(whichMod), nbins - 1);
        elseif nbins == 2,
            uncQs = median(data.(whichMod));
        end
        
        % uncertainty bins
        for u = 1:nbins,
            
            switch u
                case 1
                    trls = find(data.(whichMod) < uncQs(u));
                case nbins
                    trls = find(data.(whichMod) > uncQs(u-1));
                otherwise
                    trls = find(data.(whichMod) > uncQs(u-1) & data.(whichMod) < uncQs(u));
            end
            
            % only use correct trials!
            % trls = intersect(trls, find(data.correct == 1));
            
            % with this selection, take the trials after that
            laggedtrls = trls+lag;
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
            
            [b] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                'binomial','link','logit');
            
            % save betas
            grandavg.logistic(sj, u, lag, :) = b;
            
        end % uncertainty bin
    end % sj
end

% ========================================================= %
% normalize by their overall slope
% ========================================================= %

% grandavg.logistic = bsxfun(@minus, grandavg.logistic, grandavg.overallLogistic(:, 2));

% ========================================================= %
% combine the lag groups
% ========================================================= %

grandavg.logisticSlope = squeeze(mean(grandavg.logistic(:, :, lagGroups, 2), 3));
colors = linspecer(9);
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'pupil'));

switch grouping
    
    case 'all'
        theseSj = 1:27;
        thiscolor = colors(2, :);
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
        thiscolor = colors(8, :);
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
        thiscolor = colors(9, :);
end

stimx2 = 1:u;
hold on;

errorbar(stimx2, ...
    nanmean(grandavg.logisticSlope(theseSj, :)), ...
    nanstd(grandavg.logisticSlope(theseSj, :)) ./ sqrt(length(theseSj)), ...
    'o-', 'color', thiscolor, 'markeredgecolor', 'w', 'markerfacecolor', thiscolor);

xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins, ...
    'xticklabel', {'low', 'med', 'high'}, 'xticklabelrotation', 0);
box off;
axis tight; axis square;

cnt = 0; clear x s f
for sj = 1:length(theseSj),
    for u = 1:nbins,
        cnt = cnt + 1;
        x(cnt) = grandavg.logisticSlope(theseSj(sj), u);
        s(cnt) = sj;
        f{1}(cnt) = u; % factor of previous uncertainty
    end
end

clear stats
stats = rm_anova(x, s, f);

ymax = max( nanmean(grandavg.logisticSlope(theseSj, :)) + ...
    2* nanstd(grandavg.logisticSlope(theseSj, :)) ./ sqrt(length(theseSj)));
mysigstar([stimx2(1) stimx2(end)], [ymax ymax], stats.f1.pvalue);

switch whichMod,
    case 'baseline_pupil';
        xlabel(sprintf('Baseline pupil', 0));
    case 'decision_pupil'
        xlabel(sprintf('Pupil response'));
    case 'rt'
        xlabel(sprintf('Reaction time'));
    case 'feedback_pupil',
        xlabel('Feedback pupil');
        
        switch whichmodulator
            case 'fb-decpupil'
                xlabel('Feedback-dec pupil');
        end
    otherwise
        %  xlabel(sprintf('%s_{t-%d}', whichModTitle, whichLag));
end
ylabel({'Next trial slope'});
xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins+1, ...
    'xticklabel', {'low', 'med', 'high'}, 'xticklabelrotation', 0);
box off;

ylims = get(gca, 'ylim');
set(gca, 'ylim', [ylims(1) - 0.2*range(ylims) ylims(2)+0.2*range(ylims)]);

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4d_psychFuncShift_slope.pdf'));

end
