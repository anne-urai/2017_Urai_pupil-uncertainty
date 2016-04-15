function [] = ReactionTimeUncertainty()
% plot RT as a function of |motionstrength|, to show that the pattern
% follows model predictions of uncertainty

global mypath;
nbins = 6;

% can try this also with all subjects
subjects = 1:27;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grandavg.rt.data = nan(length(subjects), 2, nbins);

for sj = subjects,
    
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % normalization etc
    data.motionstrength = (abs(data.motionstrength));
    
    % project pupil out of RT so that we are sure these effects are independent
    data.cleanRt = projectout(data.rt, zscore(data.decision_pupil));
    assert(isreal(data.cleanRt));
    
    % loop over error and correct
    cors = [0 1];
    for corr = cors,
        
        % FIT BETAS ON THE FULL MODEL, NOT BINNED
        trls = find(data.correct == corr);
        
        % USE EVIDENCE STRENGTH TO PREDICT log(RT)
        % add pupil into the design matrix
        mdl = fitlm([data.motionstrength(trls) zscore(data.decision_pupil(trls))],  ...
            zscore(log(data.rt(trls) + 0.1)), 'robustopts', 'on');
        grandavg.rt.regline(find(sj==subjects), find(corr==cors), :) = ...
            mdl.Coefficients.Estimate;
        
        % RATHER THAN DISCRETE CATEGORIES, BIN BY RPE
        clear trls;
        trls = find(data.subjnr==sj & data.correct==corr);
        
        % get the mean pupil dilation out
        [grandavg.rpeMean(find(sj==subjects), find(corr==cors), :), ...
            grandavg.rt.data(find(sj==subjects), find(corr==cors), :), ...
            grandavg.rpeStd(find(sj==subjects), find(corr==cors), :), ...
            grandavg.rt.wgt(find(sj==subjects), find(corr==cors), :)] = ...
            divideintobins(data.motionstrength(trls), data.rt(trls), nbins);
    end
end

% PLOT
% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 8);
cols = colors([1 3], :);

%% for erlich, plot with all SJ

figure; subplot(3,3,1);
hold on;
for co = 1:2,
    plot(squeeze(grandavg.rpeMean(:, co, :))', squeeze(grandavg.rt.data(:, co, :))', ...
        '.-', 'color', cols(co, :));
end

xlabel('Evidence');
ylabel('Reaction time (s)');
ylim([0 1.2]); set(gca, 'ytick', 0:0.5:1.5);
xlim([-0.2 5.5]); set(gca, 'xtick', 0:2.75:5.5, 'xticklabel', {'weak', 'medium', 'strong' });
axis square;
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
saveas(gcf, sprintf('~/Dropbox/Meetings/RT.eps'), 'epsc');

figure; subplot(3,3,1);
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

xlabel('Evidence');
ylabel('Reaction time (s)');
ylim([0.25 0.7]); set(gca, 'ytick', [0.3 0.4 0.5 0.6]);
xlim([-0.2 5.5]); set(gca, 'xtick', 0:2.75:5.5, 'xticklabel', {'weak', 'medium', 'strong' });
axis square;
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
saveas(gcf, sprintf('~/Dropbox/Meetings/GA_RT.eps'), 'epsc');

figure; subplot(3,3,1);
% normalize each sj to mean!
for sj = subjects,
    meanRT = grandavg.rt.data(sj, :, :);
    grandavg.rt.data(sj, :, :) = grandavg.rt.data(sj, :, :) - mean(meanRT(:)) + 0.5;
end

hold on;
for co = 1:2,
    plot(squeeze(grandavg.rpeMean(:, co, :))', squeeze(grandavg.rt.data(:, co, :))', ...
        '.-', 'color', cols(co, :));
end

xlabel('Evidence');
ylabel('Normalized RT (s)');
ylim([0 1.2]); set(gca, 'ytick', 0:0.5:1.5);
xlim([-0.2 5.5]); set(gca, 'xtick', 0:2.75:5.5, 'xticklabel', {'weak', 'medium', 'strong' });
axis square;
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
saveas(gcf, sprintf('~/Dropbox/Meetings/NormalizedRT.eps'), 'epsc');


%% make the subplot next to it show the significance of the intercepts
% and slopes
close all;

% slopes
slopes         = [grandavg.rt.regline(:, 1, 2) grandavg.rt.regline(:, 2, 2)];
[~, pvalE_interc, ~, stat] = ttest(slopes(:, 1), 0, 'tail', 'both');
bf10 = t1smpbf(stat.tstat,27);
[~, pvalC_interc, ~, stat] = ttest(slopes(:, 2), 0, 'tail', 'both');
bf10 = t1smpbf(stat.tstat,27);
[~, pvalD_interc, ~, stat] = ttest(slopes(:,1), slopes(:,2));
bf10 = t1smpbf(stat.tstat,27);

% slopes
subplot(3,3,1);
hold on;
bar(1, mean(slopes(:,1)), 'FaceColor',  cols(1, :), 'EdgeColor', 'none', 'BarWidth', 0.8);
bar(2, mean(slopes(:,2)), 'FaceColor', cols(2, :), 'EdgeColor', 'none', 'BarWidth', 0.8);
h = ploterr(1:2, mean(slopes), [], std(slopes)/ sqrt(length(subjects)), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none');
axis square;
set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'}, 'ytick', [-1:0.1:1]);
ylabel('Beta weight (a.u.)');
plot(slopes', '.k-', 'linewidth', 0.2, 'markersize', 2);

mysigstar([1 2], [0.4], pvalD_interc); % difference
mysigstar(1, 0.35, pvalE_interc);
mysigstar(2, -0.35, pvalC_interc);

set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'ytick', [-0.3:0.3:0.3]);
ylim([-0.4 0.4]); xlim([0.5 2.5]); 
saveas(gcf, sprintf('~/Dropbox/Meetings/RTbetas.eps'), 'epsc');

end
