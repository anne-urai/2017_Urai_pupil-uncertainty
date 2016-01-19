function fig4_psychFuncShift_RT(lagGroups, whichmodulator, grouping, correctness)
% after talking to Redmond, split not only by pupil and choice but also
% next choice - and calculate RT!

if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('grouping', 'var'); grouping = 'all'; end
if ~exist('correctness', 'var'); correctness = []; end
warning('error', 'stats:glmfit:PerfectSeparation');

whichMod = 'decision_pupil';
nbins = 4;
subjects = 1:27;
clear grandavg;

lag = 1;
for sj = unique(subjects),
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
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
            grandavg.rt(sj, r, u, 1) = mean(thisdat.rt(thisdat.resp == resps(r)));
            grandavg.rt(sj, r, u, 2) = mean(thisdat.rt(thisdat.resp ~= resps(r)));
            
        end % uncertainty bin
    end % resp
    
    % also compute their overall RT
    %grandavg.meanRT(sj) = mean(data.rt);
    %grandavg.rt(sj, :, :, :) = grandavg.rt(sj, :, :, :) - mean(data.rt);
    
end % sj

%grandavg.rt = grandavg.rt + mean(grandavg.meanRT); % normalize
grandavg.rt = grandavg.rt * 1000; % ms

% ========================================================= %
% compute RT advantage for repeating
% ========================================================= %

% rt advantage when repeating the choice
grandavg.rtAdvantage = grandavg.rt(:, :, :, 2) - grandavg.rt(:, :, :, 1);
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
    (nanstd(grandavg.rtAdvantage(theseSj, :))) ./sqrt(length(subjects)), ...
     '.-', 'color', thiscolor, 'markerfacecolor', thiscolor, 'markersize', 12);

% low vs high
% one sample t-test with the prediction of less switching
[~, pval] = ttest(grandavg.rtAdvantage(theseSj, 1), grandavg.rtAdvantage(theseSj, end), 'tail', 'right');
[pval] = permtest(grandavg.rtAdvantage(theseSj, 1), grandavg.rtAdvantage(theseSj, end));

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
xlabel('Current trial pupil');
ylabel({'RT repetition advantage';'on next trial'});

end
