function fig4_psychFuncShift_RT(whichmodulator, grouping, correctness, nbins)
% after talking to Redmond, split not only by pupil and choice but also
% next choice - and calculate RT!

if ~exist('whichmodulator', 'var'); whichmodulator = 'rt'; end
if ~exist('grouping', 'var'); grouping = 'all'; end
if ~exist('correctness', 'var'); correctness = []; end
if ~exist('nbins', 'var'); nbins = 3; end
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
    case 'evidence'
        whichMod = 'evidence';
end

subjects = 1:27;
clear grandavg;

lag = 1;
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
        case 'evidence'
            data.evidence = abs(data.motionstrength);
    end
    
    % outcome vector need to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    
    % previous response
    resps = [0 1];
    for r = 1:2,
        
        % split into pupil bins
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
            thisdat = data(laggedtrls, :);
            
            % compute RT separately for the repeating and the alternating choice
            grandavg.rt(sj, r, u, 1) = median(thisdat.rt(thisdat.resp == resps(r)));
            grandavg.rt(sj, r, u, 2) = median(thisdat.rt(thisdat.resp ~= resps(r)));
            
        end % uncertainty bin
    end % resp
end % sj

%grandavg.rt = grandavg.rt + mean(grandavg.meanRT); % normalize
grandavg.rt = grandavg.rt * 1000; % ms

% ========================================================= %
% compute RT advantage for repeating
% ========================================================= %
% 
% hold on;
% colors = linspecer(10);
% rtswitch = grandavg.rt(:, :, :, 2);
% rtrep = grandavg.rt(:, :, :, 1);
% for u = 1:nbins,
%     for r = 1:2,
%         plot(grandavg.rt(:, r, u, 1), grandavg.rt(:, r, u, 2), '.', 'color', colors(u, :));
%     end
% end
% 
% hold on; plot([0:1:1000], 0:1000, 'k');
% xlim([0 1000]); ylim([0 1000]); axis square;

% rt advantage when repeating the choice
grandavg.rtAdvantage = grandavg.rt(:, :, :, 2) - grandavg.rt(:, :, :, 1); % faster when repeating
grandavg.rtAdvantage = grandavg.rtAdvantage(:, 1, :) + grandavg.rtAdvantage(:, 2, :) ./ 2;
grandavg.rtAdvantage = squeeze(grandavg.rtAdvantage);

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

errorbar(stimx2, (nanmean(grandavg.rtAdvantage(theseSj, :))), ...
    (nanstd(grandavg.rtAdvantage(theseSj, :))) ./sqrt(length(theseSj)), ...
    '.-', 'color', thiscolor, 'markerfacecolor', thiscolor, 'markersize', 12);

% low vs high
% one sample t-test with the prediction of less switching
[~, pval] = ttest(grandavg.rtAdvantage(theseSj, 1), grandavg.rtAdvantage(theseSj, end))
% [pval] = permtest(grandavg.rtAdvantage(theseSj, 1), grandavg.rtAdvantage(theseSj, end));

ymax = max( nanmean(grandavg.rtAdvantage(theseSj, :)) + ...
    2* nanstd(grandavg.rtAdvantage(theseSj, :)) ./ sqrt(length(theseSj)));
mysigstar([1 nbins], [ymax ymax], pval, thiscolor, 'down');

box off;
xlim([0.5 nbins+0.5]); set(gca, 'xtick', 1:nbins, ...
    'xticklabel', {'low', 'med', 'high'});
switch grouping
    case 'all'
        % ylim([0.49 0.54]);
    case 'repeat'
        %  ylim([0.5 0.55]);
        title(tit, 'color', titcolor);
    case 'switch'
        % ylim([0.45 0.53]);
        title(tit, 'color', titcolor);
end
%set(gca, 'ytick', 0.38:0.02:0.6);
axis square;
switch whichmodulator
    case 'pupil'
        xlabel('Current trial pupil');
    case 'evidence'
        xlabel('Current trial evidence');
    case 'rt'
        xlabel('Current trial RT');
end
ylabel({'RT repetition advantage';'on next trial'});

end
