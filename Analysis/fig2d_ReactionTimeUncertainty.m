function [] = f1e_RT_uncertainty()
% plot RT as a function of |motionstrength|, to show that the pattern
% follows model predictions of uncertainty

close all; clear; clc;

% load('~/Data/UvA_pupil/GrandAverage/reinforcementlearning_fmincon_fixbc.mat');
data = readtable('~/Data/UvA_pupil/CSV/2ifc_data_allsj.csv');

nbins               = 6; % bin in 5 to have comparable plots to the difficulty version?
cohs                = 1:nbins;
data.xval           = abs(data.motionstrength);
data.rpebin         = nan(size(data.xval));

% can try this also with all subjects
subjects = 1:27;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields      = {'rt'}; f = 1;
figure;
set(gcf, 'DefaultAxesFontSize', 8, 'DefaultAxesFontName', 'Helvetica');

    grandavg.(fields{f}).data = nan(length(subjects), 2, nbins);
    
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
            grandavg.(fields{f}).regline(find(sj==subjects), find(corr==cors), :) = ...
                mdl.Coefficients.Estimate;
            
            % RATHER THAN DISCRETE CATEGORIES, BIN BY RPE
            clear trls;
            trls = find(data.subjnr==sj & data.correct==corr);
            
            % get the mean pupil dilation out
            [grandavg.rpeMean(find(sj==subjects), find(corr==cors), :), ...
                grandavg.(fields{f}).data(find(sj==subjects), find(corr==cors), :), ...
                grandavg.rpeStd(find(sj==subjects), find(corr==cors), :), ...
                grandavg.(fields{f}).wgt(find(sj==subjects), find(corr==cors), :) ] = ...
                divideintobins(data.xval(trls), data.(fields{f})(trls), nbins);
            
        end
    end
    
    % PLOT
    subplot(4,4,1);
    
    % use nice shades of red and green
    cols = linspecer(3); cols = cols(2:3, :);
    
    hold on;
    for co = 1:2,
        if co == 1,
            h = ploterr(squeeze(nanmean(grandavg.rpeMean(:, co, :))), ...
                squeeze(nanmean(grandavg.(fields{f}).data(:, co, :))),...
                squeeze(nanstd(grandavg.rpeMean(:, co, :))) / sqrt(length(subjects)), ...
                squeeze(nanstd(grandavg.(fields{f}).data(:, co, :))) / sqrt(length(subjects)), ...
                'ks-', 'hhxy', 0.0011);
            set(h(1), 'color', cols(co, :), ...
                'markerfacecolor', cols(co, :),  'markersize', 3);
        elseif co == 2,
            h = ploterr(squeeze(nanmean(grandavg.rpeMean(:, co, :))), ...
                squeeze(nanmean(grandavg.(fields{f}).data(:, co, :))),...
                squeeze(nanstd(grandavg.rpeMean(:, co, :))) / sqrt(length(subjects)), ...
                squeeze(nanstd(grandavg.(fields{f}).data(:, co, :))) / sqrt(length(subjects)), ...
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
    xlabel('Stimulus difficulty');
    ylabel('Reaction time (s)');
    ylim([0.8 1.2]); set(gca, 'ytick', [0.8 1 1.2]);
    xlim([-4 4]); set(gca, 'xtick', 0:4, 'xticklabel', {'c_1', 'c_2', 'c_3', 'c_4', 'c_5'});
    offsetAxes(gca, 0.1, 0);
    set(gca, 'xcolor', 'k', 'ycolor', 'k');

    %% make the subplot next to it show the significance of the intercepts
    % and slopes
    
    if 0,
        subplot(length(fields),length(fields),length(fields)*(f-1)+2);
        
        % slopes
        slopes         = [grandavg.(fields{f}).regline(:, 1, 2) grandavg.(fields{f}).regline(:, 2, 2)];
        [~, pvalE_interc, ~, stat] = ttest(slopes(:, 1), 0, 'tail', 'both');
        bf10 = t1smpbf(stat.tstat,27)
        [~, pvalC_interc, ~, stat] = ttest(slopes(:, 2), 0, 'tail', 'both');
        bf10 = t1smpbf(stat.tstat,27)
        [~, pvalD_interc, ~, stat] = ttest(slopes(:,1), slopes(:,2));
        bf10 = t1smpbf(stat.tstat,27)

        % slopes
        hold on;
        bar(1, mean(slopes(:,1)), 'FaceColor',  cols(1, :), 'EdgeColor', 'w', 'BarWidth', 0.4);
        bar(2, mean(slopes(:,2)), 'FaceColor', cols(2, :), 'EdgeColor', 'w', 'BarWidth', 0.4);
        errorbar(1:2, mean(slopes), std(slopes)/ sqrt(length(subjects)), 'k', 'Marker', 'none', 'LineStyle', 'none');
        xlim([0.5 2.5]); set(gca, 'tickdir', 'out', 'xtick', 1:2, 'xticklabel', ...
            [] , 'ydir', 'normal', 'xticklabelrotation', 0);
        
        if f == length(fields), set(gca, 'xticklabel', {'Error', 'Correct'}); end
        ylabel('\beta');
        sigstar({[1 2]}, pvalD_interc);
        sigstar({[1,1], [2,2]}, [pvalE_interc pvalC_interc]);
        ylims = get(gca, 'ylim');  ylim([-max(abs(ylims)) max(abs(ylims))]);
    end
    

set(gca, 'xcolor', 'k', 'ycolor', 'k');

print(gcf, '-dpdf', sprintf('~/Dropbox/Writing/AnneUrai_figures/learning/Fig1e_RTuncertainty.pdf'));
end
