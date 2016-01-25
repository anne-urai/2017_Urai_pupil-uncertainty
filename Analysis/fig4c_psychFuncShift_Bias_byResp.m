function fig4c_psychFuncShift_Bias_byResp(whichmodulator, grouping, nbins)
global mypath;

if ~exist('whichmodulator', 'var'); whichmodulator = 'evidence'; end
if ~exist('grouping', 'var'); grouping = 'all'; end

% plot both the effect of pupil on overall repetition bias and show that
% this is symmetrical for both previous choices
hold on;

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
    case 'evidence'
        whichMod = 'evidence';
end

subjects = 1:27;
clear grandavg;
whichLags = 1:3; % lag 1

for lag = whichLags,
    for sj = unique(subjects),
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
            data = data(find(data.sessionnr > 1), :);

        % in this case, take out decision effects
        switch whichmodulator
            case 'fb-decpupil'
                data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
            case 'pupil-rt',
                data.decision_pupil = projectout(data.decision_pupil, data.rt);
            case 'baselinepupil'
                data.baseline_pupil = circshift(data.baseline_pupil, -1);
            case 'evidence'
                data.evidence = abs(data.motionstrength);
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
            
            % select for response and correctness
            trls = find(data.resp == resps(r));
            
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
% transform to probability saying '1'
% ========================================================= %

grandavg.logistic = exp(grandavg.logistic)./(1+exp(grandavg.logistic));

% ========================================================= %
% combine the lag groups
% ========================================================= %

grandavg.logistic = squeeze(nanmean(grandavg.logistic(:, :, :, 1, :), 4));

% ========================================================= %
% plot for all subjects
% ========================================================= %

colors = cbrewer('qual', 'Set2', 8);
colors = colors([1 3], :);
stimx2 = 1:nbins;

% split subjects based on their plain history weights
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));

hold on;
switch grouping
    
    case 'all'
        theseSj = 1:27;
    case 'repeat'
        theseSj = find(dat.response(:, 1) > 0);
    case 'switch'
        theseSj = find(dat.response(:, 1) < 0);
end

plot([1 nbins], [0.5 0.5], 'k', 'linewidth', 0.5);

for r = [1 2],
    errorbar(stimx2, ...
        squeeze(nanmean(grandavg.logistic(theseSj, r, :, 1))),   ...
        squeeze(nanstd(grandavg.logistic(theseSj, r, :, 1))) ./ sqrt(length(theseSj)), ...
        '.-', 'color', colors(r, :), 'markersize', 12);
    stimx2 = stimx2 - 0.05;
    
end

axis tight; axis square;
xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins, 'xticklabel', {'low', 'med', 'high'});
ylabel('Next trial P(choice A)');

switch whichMod
    case 'pupil'
        xlabel('Current trial pupil');
    case 'evidence'
        xlabel('Current trial evidence');
    otherwise
        xlabel(whichmodulator);
end


text(2.5, .53, 'current', 'color', colors(2, :), 'horizontalalignment', 'center');
text(2.5, 0.52, 'choice A', 'color', colors(2, :), 'horizontalalignment', 'center');

text(2.5, 0.48, 'current', 'color', colors(1, :), 'horizontalalignment', 'center');
text(2.5, 0.47, 'choice B', 'color', colors(1, :), 'horizontalalignment', 'center');

ylim([0.465 .545]);
set(gca, 'ytick', [0.47 0.5 0.53]);
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

end
