function fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, grouping, correctness)

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('grouping', 'var'); grouping = 'all'; end
if ~exist('correctness', 'var'); correctness = []; end

% plot both the effect of pupil on overall repetition bias and show that
% this is symmetrical for both previous choices
% subplot(443);

switch whichmodulator
    case 'pupil'
        whichMod = 'decision_pupil';
    case 'feedbackpupil'
        whichMod = 'feedback_pupil';
    case 'fb-decpupil'
        % see below for projecting out
        whichMod = 'feedback_pupil';
    case 'rt'
        whichMod = 'rt';
    case 'pupil-rt';
        whichMod = 'decision_pupil';
end

nbins = 3;
subjects = 1:27;
clear grandavg;
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
        end
        
        % outcome vector need to be 0 1 for logistic regression
        data.resp(data.resp == -1) = 0;
        
        % get an overall logistic fit for normalization
        grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, ...
            data.resp, ...
            'binomial','link','logit');
        
        % previous response
        resps = [0 1];
        for r = 1:2,
            
            if isempty(correctness),
                % select for response and correctness
                trls = find(data.resp == resps(r));
            else
                trls = find(data.resp == resps(r) & data.correct == correctness);
            end
            
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
            
            % get an overall logistic fit
            grandavg.respBiasNoPupilSplit(sj, r, lag, :) = glmfit(thisdat.motionstrength, thisdat.resp, ...
                'binomial','link','logit');
            
            %% split into pupil bins
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
                
                % take subset
                if ~isempty(correctness),
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
                
                [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                    'binomial','link','logit');
                
                % save betas
                grandavg.logistic(sj, r, u, lag, :) = b;
                
            end % uncertainty bin
        end % resp
    end % sj
end

% ========================================================= %
% normalize by their overall bias
% ========================================================= %

grandavg.logistic = bsxfun(@minus, grandavg.logistic, grandavg.overallLogistic(:, 1));
grandavg.respBiasNoPupilSplit = bsxfun(@minus, grandavg.respBiasNoPupilSplit, grandavg.overallLogistic(:, 1));

% ========================================================= %
% make one plot with repetition (rather than response) bias
% ========================================================= %

resp1 = squeeze(grandavg.logistic(:, 1, :, :, 1));
resp2 = squeeze(grandavg.logistic(:, 2, :, :, 1));

grandavg.logisticRep = (-resp1 + resp2) ./ 2;

resp1 = squeeze(grandavg.respBiasNoPupilSplit(:, 1, :, 1));
resp2 = squeeze(grandavg.respBiasNoPupilSplit(:, 2, :, 1));

grandavg.repetitionBias = (-resp1 + resp2) ./ 2;

% ========================================================= %
% combine the lag groups
% ========================================================= %

grandavg.logisticRep = squeeze(nanmean(grandavg.logisticRep(:, :, lagGroups), 3));
grandavg.logistic = squeeze(nanmean(grandavg.logistic(:, :, :, lagGroups, :), 4));
grandavg.respBiasNoPupilSplit = squeeze(nanmean(grandavg.respBiasNoPupilSplit(:, :, lagGroups, 1), 3));
grandavg.repetitionBias  = squeeze(nanmean(grandavg.repetitionBias(:, lagGroups), 2));

% ========================================================= %
% plot
% ========================================================= %

colors = linspecer(9);
% split subjects based on their plain history weights
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
hold on;

switch grouping
    case 'all'
        theseSj = 1:27;
        tit = 'All subjects';
        titcolor = colors(7, :);
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
        titcolor = colors(8, :);
        tit = 'Repeaters';
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
        titcolor = colors(9, :);
        tit = 'Switchers';
end

if isempty(correctness),
    thiscolor = 'k';
else
    % rather define color by correctness
    switch correctness
        case 1
            thiscolor = colors(3,:);
            plot([1 nbins], [0 0], 'k', 'linewidth', 0.5);
            
        case 0
            thiscolor = colors(1,:);
    end
end

stimx2 = 1:u;
if correctness == 0,
    stimx2 = stimx2 + 0.05;
elseif correctness == 1,
    stimx2 = stimx2 - 0.05;
end

errorbar(stimx2, ...
    nanmean(grandavg.logisticRep(theseSj, :)), ...
    nanstd(grandavg.logisticRep(theseSj, :)) ./ sqrt(length(theseSj)), ...
    '.-', 'color', thiscolor, 'markerfacecolor', thiscolor, 'markersize', 12);

axis tight; axis square;

% low vs high
[~, pval(1)] = ttest(grandavg.logisticRep(theseSj, 1), grandavg.logisticRep(theseSj, 3));
ymax = max( nanmean(grandavg.logisticRep(theseSj, :)) + ...
    2* nanstd(grandavg.logisticRep(theseSj, :)) ./ sqrt(length(theseSj)));
%mysigstar([1 3], [ymax ymax], pval(1));

if 0,
    % only test for a repetition bias when we haven't split the subjects
    switch grouping
        case 'all'
            [~, pval(3)] = ttest(grandavg.repetitionBias(theseSj));
            max( nanmean(grandavg.repetitionBias(theseSj)) + ...
                nanstd(grandavg.repetitionBias(theseSj)) ./ sqrt(length(theseSj)));
            mysigstar(nbins+1, ymax, pval(3));
    end
end

box off;

xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins, ...
    'xticklabel', []);
title(tit, 'color', titcolor);

switch grouping
    
    case 'all'
        ylim([-0.05 0.2]);
        
    case 'repeat'
        ylim([-0.05 0.4]);
        
    case 'switch'
        ylim([-0.4 0.05]);
        
        % xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins, ...
        %     'xticklabel', {'low', 'med', 'high'});
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
        end
end
ylabel('Repetition');


end
