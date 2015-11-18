function [] = f2c_PupilUncertaintyCorrelation_Feedback()
% using the timewindow that is indicated in the regression timecourse plot,
% show the shape of the pupil vs motionstrength pattern

close all; clear; clc;
RTstratification    = true; % include RT in the model, do stratification on bins
RTbinsize           = 0.01; % the larger the binsize, the more trials we can keep (in seconds)

% get all data
data = readtable('~/Data/UvA_pupil/CSV/2ifc_data_allsj.csv');
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens
subjects = 1:27; % for this analysis, use all SJ!

nbins       = 6; % bin in 5 to have comparable plots to the difficulty version?
data.xval   = abs(data.motionstrength);
data.rpebin = nan(size(data.xval));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields      = {'feedback_pupil', };
fieldnames  = {'Pupil response'};
fieldunits  = {'(% signal change)'};
figure;

for f = 1:length(fields),
    
    grandavg.(fields{f}).data = nan(length(subjects), 2, nbins);
    
    for sj = subjects,
        
        % loop over error and correct
        cors = [0 1];
        for corr = cors,
            
            % FIT BETAS ON THE FULL MODEL, NOT BINNED
            trls = find(data.subjnr == sj & data.correct == corr);
            
            if RTstratification,
                % include RT as a regressor
                mdl = fitlm([zscore(abs(data.xval(trls))) zscore(data.rt(trls))],  ...
                    zscore(data.(fields{f})(trls)));
            else
                % dont include RT as a regressor?
                mdl = fitlm(([abs(data.xval(trls))]),  ...
                    (data.(fields{f})(trls)));
            end
            
            % SAVE BETAS FOR THIS PARTICIPANT
            grandavg.(fields{f}).regline(find(sj==subjects), find(corr==cors), :) = ...
                mdl.Coefficients.Estimate;
            
            % PRE-STRATIFY THE RT DISTRIBUTIONS
            if 0,
                if RTstratification && f == 1,
                    
                    tmpdata = data;
                    % DIVIDE INTO BINS OF EXPECTED REWARD
                    tmpdata.rpebin(trls) = quantileIdx(data.xval(trls), nbins);
                    
                    for coh = cohs,
                        trls = find(tmpdata.subjnr==sj & tmpdata.correct==corr & tmpdata.rpebin==coh);
                        rts{coh} = tmpdata.rt(trls);
                    end
                    
                    % get indices we want to keep for each distribution
                    [idx_keep, idx_reject] = commonDistributions(RTbinsize, rts);
                    
                    % now select this specific subset of data
                    for coh = cohs,
                        trls = find(tmpdata.subjnr==sj & tmpdata.correct==corr & tmpdata.rpebin==coh);
                        tmpdata(trls(idx_reject{coh}), :) = [];
                        fprintf('sj %d, coh %d, corr %d, removing %d trials \n', sj, coh, corr, numel(idx_reject{coh}));
                    end
                else
                    tmpdata = data;
                    
                end
            else
                tmpdata = data;
                
            end
            
            % RATHER THAN DISCRETE CATEGORIES, BIN BY motionenergy
            clear trls;
            trls = find(tmpdata.subjnr==sj & tmpdata.correct==corr);
            
            % get the mean pupil dilation out
            [grandavg.xMean(find(sj==subjects), find(corr==cors), :), ...
                grandavg.(fields{f}).data(find(sj==subjects), find(corr==cors), :), ...
                grandavg.xStd(find(sj==subjects), find(corr==cors), :), ...
                grandavg.(fields{f}).wgt(find(sj==subjects), find(corr==cors), :)] = ...
                divideintobins(tmpdata.xval(trls), tmpdata.(fields{f})(trls), nbins);
            
        end
    end
    
    % PLOT
    subplot(4,4,1);
    % use nice shades of red and green
    cols = linspecer(3); cols = cols(2:3, :);
    
    % slopes
    slopes       = [grandavg.(fields{f}).regline(:, 1, 2) grandavg.(fields{f}).regline(:, 2, 2)];
    [~, pval(1), ~, stat] = ttest(slopes(:, 1), 0, 'tail', 'both');
    bf10 = t1smpbf(stat.tstat,27)
    [~, pval(2), ~, stat] = ttest(slopes(:, 2), 0, 'tail', 'both');
    bf10 = t1smpbf(stat.tstat,27)
    [~, pval(3), ~, stat] = ttest(slopes(:,1), slopes(:,2));
    bf10 = t1smpbf(stat.tstat,27)
    
    hold on;
    for co = 1:2,
        
        % use double error bars
        if co == 1,
            h = ploterr(squeeze(mean(grandavg.xMean(:, co, :))), ...
                squeeze(nanmean(grandavg.(fields{f}).data(:, co, :))), ...
                squeeze(nanstd(grandavg.xMean(:, co, :))) / sqrt(length(subjects)), ...
                squeeze(nanstd(grandavg.(fields{f}).data(:, co, :))) / sqrt(length(subjects)), ...
                'ks-',  'hhxy', 0.001);
            set(h(1), 'color', cols(co, :), ...
                'markerfacecolor', cols(co, :),  'markersize', 3);
        elseif co == 2,
            h = ploterr(squeeze(mean(grandavg.xMean(:, co, :))), ...
                squeeze(nanmean(grandavg.(fields{f}).data(:, co, :))), ...
                squeeze(nanstd(grandavg.xMean(:, co, :))) / sqrt(length(subjects)), ...
                squeeze(nanstd(grandavg.(fields{f}).data(:, co, :))) / sqrt(length(subjects)), ...
                'ko-',  'hhxy', 0.001);
            
            set(h(1), 'color', cols(co, :), ...
                'markerfacecolor', cols(co, :),  'markersize', 12, 'marker', '.');
        end
        set(h(2), 'color', cols(co, :) + 0.05);
        set(h(3), 'color', cols(co, :) + 0.05);
    end
    
    if f == length(fields),
        xlabel('|\Delta| motion strength (a.u.)');
    else
        set(gca, 'xticklabel', []);
    end
    set(gca, 'xtick', 0:4, 'xticklabel',  {'c_1', 'c_2', 'c_3', 'c_4', 'c_5'});
    
  set(gca, 'ytick', [-.5:0.5:2]);
    ylabel({fieldnames{f}; fieldunits{f}});
   offsetAxes(gca, 0.1, 0); 
    set(gca, 'xcolor', 'k', 'ycolor', 'k');

    % make a barplot
    % subplot(length(fields),length(fields), length(fields)*(f-1)+2);
    subplot(4,5,3);
    h = ploterr(1:2, squeeze(mean(slopes)), [], squeeze(std(slopes)) / sqrt(length(subjects)), 'k.', 'hhxy', 0.00001);
    set(h(1), 'marker', 'none');
    
    hold on;
    bar(1, mean(slopes(:, 1)), 'facecolor', cols(1, :), 'edgecolor', 'none', 'barwidth', 0.4);
    bar(2, mean(slopes(:, 2)), 'facecolor', cols(2, :), 'edgecolor', 'none', 'barwidth', 0.4);
    s1 = sigstar({[1 1]}, pval(1), 0); %set(s1(2), 'position', [1 0.01 0]);
    s1 = sigstar({[2 2]}, pval(2), 0); %set(s1(2), 'position', [2 -0.01 0]);
    s1 = sigstar({[1 2]}, pval(3), 0); 
    
    set(gca, 'xtick', [0.5 1 2 2.5]);
    try ylims = get(gca, 'ylim');
        set(gca, 'ytick', [-.05 -0.025 0 .025]);     end
    
    offsetAxes; 
   box off;
    if f == length(fields),
        set(gca, 'xtick', 1:2, 'xticklabel', {'Error', 'Correct'});
    else
        set(gca, 'xtick', 1:2, 'xticklabel', []);
    end
    % add sigstars
    ylabel('Beta weights (a.u.)');
    
end
set(gca, 'xcolor', 'k', 'ycolor', 'k');
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/learning/Fig3a_scalarPupilFeedback.pdf'));

savefast('~/Data/UvA_pupil/GrandAverage/grandavg_pupil_uncertainty_feedback.mat', 'grandavg');

end
