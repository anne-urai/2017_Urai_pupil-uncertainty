function fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, grouping, correctness)
% plot both the effect of pupil on overall repetition bias and show that
% this is symmetrical for both previous choices

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('grouping', 'var'); grouping = 'all'; end
if ~exist('correctness', 'var'); correctness = []; end

%%
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
    case 'pupil-rt'
        whichMod = 'decision_pupil';
    case 'baselinepupil'
        whichMod = 'baseline_pupil';
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
            case 'pupil-rt',
                data.decision_pupil = projectout(data.decision_pupil, data.rt);
            case 'baselinepupil'
                data.baseline_pupil = circshift(data.baseline_pupil, -1);
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
            
            if ~isempty(correctness),
                % only use correct trials!
                trls = intersect(trls, find(data.correct == correctness));
            end
            
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
colors = cbrewer('qual', 'Set1', 9);

% split subgroups by their plain history effects
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));

switch grouping
    case 'all'
        theseSj = 1:27;
        tit = 'All subjects';
        titcolor = colors(7, :);
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
        titcolor = colors(2, :);
        tit = 'Repeaters';
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
        titcolor = colors(5, :);
        tit = 'Switchers';
end

stimx2 = 1:u;
if correctness == 0,
    stimx2 = stimx2 + 0.05;
elseif correctness == 1,
    stimx2 = stimx2 - 0.05;
end
hold on;

if isempty(correctness),
    thiscolor = colors(9, :);
else
    % rather define color by correctness
    switch correctness
        case 1
            thiscolor = colors(3,:);
        case 0
            thiscolor = colors(1,:);
    end
end

errorbar(stimx2, ...
    nanmean(grandavg.logisticSlope(theseSj, :)), ...
    nanstd(grandavg.logisticSlope(theseSj, :)) ./ sqrt(length(theseSj)), ...
    '.-', 'color', thiscolor, 'markersize', 12);

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

box off;
% title(tit, 'color', titcolor);

switch grouping
    
    case 'all'
        ylim([-0.05 0.2]);
        ylabel('Next trial slope');
        
    case 'repeat'
        ylim([-0.05 0.4]);
        
    case 'switch'
        ylim([-0.4 0.05]);
        
        xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins, ...
            'xticklabel', {'low', 'med', 'high'});
end
axis tight;
xlim([0.5 nbins+0.5]);
ylim([0.65 0.85]); set(gca, 'ytick', [0.6 0.7 0.8]);
xlabel('Current trial pupil');

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4d_psychFuncShift_slope.pdf'));

end
