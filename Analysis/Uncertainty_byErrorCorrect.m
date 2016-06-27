function [b] = Uncertainty_byErrorCorrect(field)
% using the timewindow that is indicated in the regression timecourse plot,
% show the shape of the pupil vs motionstrength pattern

global mypath;

% get all data
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens
subjects = 1:27; % for this analysis, use all SJ!

nbins       = 6; % bin in 5 to have comparable plots to the difficulty version?
data.xval   = abs(data.motionstrength);
data.rpebin = nan(size(data.xval));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grandavg.(field).data = nan(length(subjects), 2, nbins);

for sj = subjects,
    
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % normalization etc
    data.motionstrength = (abs(data.motionstrength));
    
    % loop over error and correct
    cors = [0 1];
    for corr = cors,
        
        % RATHER THAN DISCRETE CATEGORIES, BIN BY motionenergy
        clear trls;
        trls = find(data.subjnr==sj & data.correct==corr);
        
        % get the mean pupil dilation out
        [grandavg.xMean(find(sj==subjects), find(corr==cors), :), ...
            grandavg.(field).data(find(sj==subjects), find(corr==cors), :), ...
            grandavg.xStd(find(sj==subjects), find(corr==cors), :), ...
            grandavg.(field).wgt(find(sj==subjects), find(corr==cors), :)] = ...
            divideintobins(data.motionstrength(trls), data.(field)(trls), nbins);
    end
    
    % project RT out of the pupil and vice versa
    switch field
        case 'rt'
            data.(field) = projectout(data.rt, zscore(data.decision_pupil));
        case 'decision_pupil'
            data.(field) = projectout(data.decision_pupil, zscore(log(data.rt + 0.1)));
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
        squeeze(nanmean(grandavg.(field).data(:, co, :))), ...
        squeeze(nanstd(grandavg.xMean(:, co, :))) / sqrt(length(subjects)), ...
        squeeze(nanstd(grandavg.(field).data(:, co, :))) / sqrt(length(subjects)), ...
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
xlim([-0.2 5.6]); set(gca, 'xtick', 0:2.75:5.5, 'xticklabel',  {'weak', 'medium', 'strong'}, 'xminortick', 'off');

switch field
    case 'decision_pupil'
        ylim([0.2 0.61]); set(gca, 'ytick', [0.2 0.4 0.6]);
        ylabel('Pupil response (z)');
    case 'rt'
        ylim([0.28 0.6]); set(gca, 'ytick', [0.3 0.4 0.5 0.6]);
        ylabel('Reaction time (s)');
end

legend([handles{:}], {'error', 'correct'}, 'location', 'southwest'); legend boxoff;
set(gca, 'xcolor', 'k', 'ycolor', 'k');

savefast(sprintf('%s/Data/GrandAverage/grandavg_pupil_uncertainty.mat', mypath), 'grandavg');

% output
b = grandavg.(field).regline;

end
