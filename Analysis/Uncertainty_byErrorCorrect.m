function [b, bint] = Uncertainty_byErrorCorrect(field, nbins)
% using the timewindow that is indicated in the regression timecourse plot,
% show the shape of the pupil vs motionstrength pattern

global mypath;

% get all data
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens
subjects = 1:27; % for this analysis, use all SJ!

if ~exist('nbins', 'var'); nbins = 6; end
data.xval   = abs(data.motionstrength);
data.rpebin = nan(size(data.xval));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grandavg.(field).data = nan(length(subjects), 2, nbins);

for sj = subjects,
    
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % if looking at the baseline, get the one from the next trial
    switch field
        case 'baseline_pupil'
            data.baseline_pupil = circshift(zscore(data.baseline_pupil), -1);
    end
    
    % normalization etc
    data.motionstrength = (abs(data.motionstrength));
    
    % loop over error and correct
    cors = [0 1];
    for corr = cors,
        
        % RATHER THAN DISCRETE CATEGORIES, BIN BY motionenergy
        clear trls;
        trls = find(data.subjnr==sj & data.correct==corr);
        
        switch field
            case 'rt'
                summaryFunc = @nanmedian;
                distFunc    = @nanstd;
            otherwise
                summaryFunc = @nanmean;
                distFunc    = @naniqr;
        end
        
        % get the mean pupil dilation out
        [grandavg.xMean(find(sj==subjects), find(corr==cors), :), ...
            grandavg.(field).data(find(sj==subjects), find(corr==cors), :), ...
            grandavg.xStd(find(sj==subjects), find(corr==cors), :), ...
            grandavg.(field).wgt(find(sj==subjects), find(corr==cors), :)] = ...
            divideintobins(data.motionstrength(trls), data.(field)(trls), nbins, [], summaryFunc);
    end
    
    % project RT out of the pupil and vice versa
    switch field
        case 'rt'
            data.(field) = projectout(zscore(data.rtNorm), zscore(data.decision_pupil));
        case 'decision_pupil'
            data.(field) = projectout(data.decision_pupil, zscore(data.rtNorm));
    end
    
    % loop over error and correct
    cors = [0 1];
    for corr = cors,
        
        % FIT BETAS ON THE FULL MODEL, NOT BINNED
        trls = find(data.subjnr == sj & data.correct == corr);
        
        % include RT as a regressor
        mdl = fitlm(zscore(data.motionstrength(trls)),  ...
            zscore(data.(field)(trls)));
        
        % SAVE BETAS FOR THIS PARTICIPANT
        grandavg.(field).regline(find(sj==subjects), find(corr==cors), :) = ...
            mdl.Coefficients.Estimate;
        bint = mdl.coefCI; coef = mdl.Coefficients.Estimate;
        grandavg.(field).bint(find(sj==subjects), find(corr==cors), :) = bint(2,2) - coef(2);
    end
    
end

% PLOT
% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 8);
cols = colors([1 2], :);

hold on;
markers = {'^', '.'}; markersizes = [4 14];
for co = 1:2,
    
    % use double error bars
    h = ploterr(squeeze(mean(grandavg.xMean(:, co, :))), ...
        squeeze(mean(grandavg.(field).data(:, co, :))), ...
        squeeze(std(grandavg.xMean(:, co, :))) / sqrt(length(subjects)), ...
        squeeze(std(grandavg.(field).data(:, co, :))) / sqrt(length(subjects)), ...
        'k-',  'abshhxy', 0);
    
    set(h(1), 'color', cols(co, :), ...
        'markersize', markersizes(co), ...
        'marker', markers{co});
    if co == 1,
        set(h(1), 'markerfacecolor', 'w', 'markeredgecolor', cols(co, :));
    end
    set(h(2), 'color', cols(co, :));
    set(h(3), 'color', cols(co, :));
    handles{co} = h(1);
end

xlabel('Evidence');
axis square;
xlim([-0.2 5.6]); set(gca, 'xtick', 0:2.75:5.5, ...
    'xticklabel',  {'weak', 'medium', 'strong'}, 'xminortick', 'off');

switch field
    case 'decision_pupil'
        ylim([0.2 0.61]); set(gca, 'ytick', [0.2 0.4 0.6]);
        ylabel('Pupil response (z)');
        savefast(sprintf('%s/Data/GrandAverage/grandavg_pupil_uncertainty.mat', mypath), 'grandavg');
    case 'baseline_pupil'
        ylim([-0.1 0.2]); set(gca, 'ytick', -1:0.1:1);
        ylabel('Next trial baseline pupil (z)');
    case 'rt'
        ylim([0.23 0.56]); set(gca, 'ytick', 0.25:0.15:0.55);
        ylabel('Reaction time (s)');
end

legend([handles{:}], {'error', 'correct'}, 'location', 'southwest'); legend boxoff;
set(gca, 'xcolor', 'k', 'ycolor', 'k');

% output
b = grandavg.(field).regline(:, :, 2);
bint = grandavg.(field).bint;
end
