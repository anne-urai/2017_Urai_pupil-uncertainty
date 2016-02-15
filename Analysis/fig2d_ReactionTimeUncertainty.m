function [] = fig2d_ReactionTimeUncertainty()
% plot RT as a function of |motionstrength|, to show that the pattern
% follows model predictions of uncertainty

global mypath;

data                = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
nbins               = 6;
data.xval           = abs(data.motionstrength);
data.rpebin         = nan(size(data.xval)); % preallocate

% can try this also with all subjects
subjects = 1:27;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grandavg.rt.data = nan(length(subjects), 2, nbins);

for sj = subjects,
    
    % loop over error and correct
    cors = [0 1];
    for corr = cors,
        
        % FIT BETAS ON THE FULL MODEL, NOT BINNED
        trls = find(data.subjnr == sj & data.correct == corr);
        
        % USE EVIDENCE STRENGTH TO PREDICT log(RT)
        mdl = fitlm(zscore(abs(data.motionstrength(trls))),  ...
            zscore(log(data.rt(trls)+0.1)));
        grandavg.rt.regline(find(sj==subjects), find(corr==cors), :) = ...
            mdl.Coefficients.Estimate;
        
        % RATHER THAN DISCRETE CATEGORIES, BIN BY RPE
        clear trls;
        trls = find(data.subjnr==sj & data.correct==corr);
        
        % get the mean pupil dilation out
        [grandavg.rpeMean(find(sj==subjects), find(corr==cors), :), ...
            grandavg.rt.data(find(sj==subjects), find(corr==cors), :), ...
            grandavg.rpeStd(find(sj==subjects), find(corr==cors), :), ...
            grandavg.rt.wgt(find(sj==subjects), find(corr==cors), :) ] = ...
            divideintobins(data.xval(trls), data.rt(trls), nbins);
        
    end
end

% PLOT
% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 8);
cols = colors([1 3], :);

hold on;
for co = 1:2,
    h = ploterr(squeeze(nanmean(grandavg.rpeMean(:, co, :))), ...
        squeeze(nanmean(grandavg.rt.data(:, co, :))),...
        squeeze(nanstd(grandavg.rpeMean(:, co, :))) / sqrt(length(subjects)), ...
        squeeze(nanstd(grandavg.rt.data(:, co, :))) / sqrt(length(subjects)), ...
        'k-', 'abshhxy', 0);
    set(h(1), 'color', cols(co, :), 'markersize', 12, 'marker', '.');
    set(h(1), 'color', cols(co, :));
    set(h(2), 'color', cols(co, :) + 0.05);
    set(h(3), 'color', cols(co, :) + 0.05);
end

% set(gca, 'box', 'off', 'tickdir', 'out', 'xtick', cohs);
xlabel('Evidence');
ylabel('Reaction time (s)');
ylim([0.25 0.61]); set(gca, 'ytick', [0.3 0.4 0.5 0.6]);
xlim([-0.2 5.5]); set(gca, 'xtick', 0:2.75:5.5, 'xticklabel', {'weak', 'medium', 'strong' });
axis square;
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

%% make the subplot next to it show the significance of the intercepts
% and slopes

subplot(4,6,23);
% slopes
slopes         = [grandavg.rt.regline(:, 1, 2) grandavg.rt.regline(:, 2, 2)];
[~, pvalE_interc, ~, stat] = ttest(slopes(:, 1), 0, 'tail', 'both');
bf10 = t1smpbf(stat.tstat,27);
[~, pvalC_interc, ~, stat] = ttest(slopes(:, 2), 0, 'tail', 'both');
bf10 = t1smpbf(stat.tstat,27);
[~, pvalD_interc, ~, stat] = ttest(slopes(:,1), slopes(:,2));
bf10 = t1smpbf(stat.tstat,27);

% slopes
hold on;
bar(1, mean(slopes(:,1)), 'FaceColor',  cols(1, :), 'EdgeColor', 'none', 'BarWidth', 0.8);
bar(2, mean(slopes(:,2)), 'FaceColor', cols(2, :), 'EdgeColor', 'none', 'BarWidth', 0.8);
h = ploterr(1:2, mean(slopes), [], std(slopes)/ sqrt(length(subjects)), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none');
axis square;
set(gca, 'xtick', [1 2], 'xticklabel', {'error', 'correct'}, 'ytick', [-1:0.1:1]);
ylabel('Beta evidence');
mysigstar([1 2], [0.1 0.1], pvalD_interc);
mysigstar([1], [0.01 0.01], pvalE_interc);
mysigstar([2], [-0.1 -0.1], pvalC_interc);

set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
ylim([-0.3 0.1]); xlim([0.1 4.5]); 

end
