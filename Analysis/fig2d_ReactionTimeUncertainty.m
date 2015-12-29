function [] = fig2d_ReactionTimeUncertainty()
% plot RT as a function of |motionstrength|, to show that the pattern
% follows model predictions of uncertainty

% load('~/Data/pupilUncertainty/GrandAverage/reinforcementlearning_fmincon_fixbc.mat');
data = readtable('~/Data/pupilUncertainty/CSV/2ifc_data_allsj.csv');

nbins               = 6; % bin in 5 to have comparable plots to the difficulty version?
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
        
        % dont include RT as a regressor?
        mdl = fitlm(([abs(data.motionstrength(trls))]),  ...
            (data.rt(trls)));
        
        % SAVE BETAS FOR THIS PARTICIPANT
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
cols = linspecer(3); cols = cols(2:3, :);

hold on;
for co = 1:2,
    if co == 1,
        h = ploterr(squeeze(nanmean(grandavg.rpeMean(:, co, :))), ...
            squeeze(nanmean(grandavg.rt.data(:, co, :))),...
            squeeze(nanstd(grandavg.rpeMean(:, co, :))) / sqrt(length(subjects)), ...
            squeeze(nanstd(grandavg.rt.data(:, co, :))) / sqrt(length(subjects)), ...
            'ks-', 'hhxy', 0.0011);
        set(h(1), 'color', cols(co, :), ...
            'markerfacecolor', cols(co, :),  'markersize', 3);
    elseif co == 2,
        h = ploterr(squeeze(nanmean(grandavg.rpeMean(:, co, :))), ...
            squeeze(nanmean(grandavg.rt.data(:, co, :))),...
            squeeze(nanstd(grandavg.rpeMean(:, co, :))) / sqrt(length(subjects)), ...
            squeeze(nanstd(grandavg.rt.data(:, co, :))) / sqrt(length(subjects)), ...
            'k-', 'hhxy', 0.0011);
        set(h(1), 'color', cols(co, :), ...
            'markerfacecolor', cols(co, :),  'markersize', 12, 'marker', '.');
    end
    
    set(h(1), 'color', cols(co, :), ...
        'markerfacecolor', cols(co, :));
    set(h(2), 'color', cols(co, :) + 0.05);
    set(h(3), 'color', cols(co, :) + 0.05);
end

% set(gca, 'box', 'off', 'tickdir', 'out', 'xtick', cohs);
xlabel('Task difficulty');
ylabel('Reaction time (s)');
ylim([0.3 0.7]); set(gca, 'ytick', [0.3 0.5 0.7]);
xlim([0 5.6]); set(gca, 'xtick', 0:2.75:5.5, 'xticklabel', {'hard', 'medium', 'easy' });
axis square;
set(gca, 'xcolor', 'k', 'ycolor', 'k');
%legend('Error','Correct'); legend boxoff;

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
bar(1, mean(slopes(:,1)), 'FaceColor',  cols(1, :), 'EdgeColor', 'w', 'BarWidth', 0.4);
bar(2, mean(slopes(:,2)), 'FaceColor', cols(2, :), 'EdgeColor', 'w', 'BarWidth', 0.4);
errorbar(1:2, mean(slopes), std(slopes)/ sqrt(length(subjects)), 'k', 'Marker', 'none', 'LineStyle', 'none');
xlim([0.5 2.5]); set(gca, 'tickdir', 'out', 'xtick', 1:2, 'xticklabel', ...
    [] , 'ydir', 'normal', 'xticklabelrotation', 0);

axis tight; 
set(gca, 'xticklabel', {'Error', 'Correct'}); 
ylabel('Beta task difficulty');
sigstar({[1 2]}, pvalD_interc);
sigstar({[1,1], [2,2]}, [pvalE_interc pvalC_interc]);
set(gca, 'xcolor', 'k', 'ycolor', 'k');
ylim([-0.05 0.08]);

end
