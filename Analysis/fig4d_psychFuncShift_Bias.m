function fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, grouping, correctness)

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'baseline'; end
if ~exist('grouping', 'var'); grouping = 'all'; end
if ~exist('correctness', 'var'); correctness = []; end
warning('error', 'stats:glmfit:PerfectSeparation');

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
    case 'baselinepupil'
        whichMod = 'baseline_pupil';
end

nbins = 4;
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
            case 'baselinepupil'
                data.baseline_pupil = circshift(data.baseline_pupil, -1);
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
            
            % fit an overall bias term
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
                
                try
                    [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                        'binomial','link','logit');
                catch
                    warning('putting zero in betas!');
                    b = nan(size(b));
                end
                
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

% ========================================================= %
% transform to probability saying '1'
% ========================================================= %

grandavg.logistic = exp(grandavg.logistic)./(1+exp(grandavg.logistic));

% ========================================================= %
% make one plot with repetition (rather than response) bias
% ========================================================= %

resp1 = -1 * (squeeze(grandavg.logistic(:, 1, :, :, 1)) - 0.5) + 0.5;
resp2 = squeeze(grandavg.logistic(:, 2, :, :, 1));

% since this is centred at 0.5, treat it that way
grandavg.logisticRep = (resp1 + resp2) ./ 2;

% ========================================================= %
% combine the lag groups
% ========================================================= %

grandavg.logisticRep = squeeze(nanmean(grandavg.logisticRep(:, :, lagGroups), 3));
grandavg.logistic = squeeze(nanmean(grandavg.logistic(:, :, :, lagGroups, :), 4));

% ========================================================= %
% plot
% ========================================================= %

colors = cbrewer('qual', 'Set1', 9);
% split subjects based on their plain history weights
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
hold on;

load('~/Data/pupilUncertainty/GrandAverage/sjcolormap.mat');
switch grouping
    case 'all'
        theseSj = 1:27;
        tit = 'All subjects';
        titcolor = 'k';
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
        titcolor = mycolmap(11,:);
        tit = 'Repeaters';
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
        titcolor = mycolmap(16,:);
        tit = 'Alternators';
end

if isempty(correctness),
    thiscolor = [0.1 0.1 0.1];
    plot([1 nbins], [0.5 0.5], 'k', 'linewidth', 0.2);
    
else
    % rather define color by correctness
    switch correctness
        case 1
            thiscolor = colors(3,:);
            plot([1 nbins], [0.5 0.5], 'k', 'linewidth', 0.2);
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

axis tight; 

% low vs high
% one sample t-test with the prediction of less switching
[~, pval] = ttest(grandavg.logisticRep(theseSj, 1), grandavg.logisticRep(theseSj, end), 'tail', 'right');
[pval] = permtest(grandavg.logisticRep(theseSj, 1), grandavg.logisticRep(theseSj, end));

%if isempty(correctness) || correctness == 0,
ymax = max( nanmean(grandavg.logisticRep(theseSj, :)) + ...
    2* nanstd(grandavg.logisticRep(theseSj, :)) ./ sqrt(length(theseSj)));
mysigstar([1 nbins], [ymax ymax], pval, thiscolor, 'down');

%elseif correctness == 1,
%    ymax = min( nanmean(grandavg.logisticRep(theseSj, :)) - ...
%        2* nanstd(grandavg.logisticRep(theseSj, :)) ./ sqrt(length(theseSj)));
%    mysigstar([1 3], [ymax ymax], pval, thiscolor, 'up');
%end

box off;
xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins, ...
    'xticklabel', {'low', 'med', 'high'});
switch grouping
    case 'all'
        ylim([0.49 0.54]);
    case 'repeat'
        ylim([0.5 0.55]);
        title(tit, 'color', titcolor);
    case 'switch'
        ylim([0.45 0.53]);
       title(tit, 'color', titcolor);
end
set(gca, 'ytick', 0.38:0.02:0.6);
axis square;
ylabel('Next trial P(repeat)');
xlabel('Current trial pupil');

end
