function fig3_BlinkSaccade
clear all; clc; close all;
addpath('~/Documents/fieldtrip');
ft_defaults;
subjects = 1:27;

load('~/Data/pupilUncertainty/GrandAverage/pupilgrandaverage.mat');
data = readtable('~/Data/pupilUncertainty/CSV/2ifc_data_allsj.csv');
data.blink = nan(size(data.decision_pupil));
data.sacc = nan(size(data.decision_pupil));

% ==================================================================
% get all data
% ==================================================================

% append all the mean timecourses per condition
for sj = unique(subjects),
    
    % get all timelock
    blinkchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'Blinks')==1);
    alltimelock.blink(sj, 1, :) = cat(2, squeeze(nanmean(pupilgrandavg.timelock{sj}(1).lock.trial(:, blinkchan, :)))', ...
        squeeze(sum(pupilgrandavg.timelock{sj}(2).lock.trial(:, blinkchan, :)))', ...
        squeeze(sum(pupilgrandavg.timelock{sj}(3).lock.trial(:, blinkchan, :)))', ...
        squeeze(sum(pupilgrandavg.timelock{sj}(4).lock.trial(:, blinkchan, :)))');
    
    % get all timelock
    saccchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'Saccades')==1);
    alltimelock.sacc(sj, 1, :) = cat(2, squeeze(sum(pupilgrandavg.timelock{sj}(1).lock.trial(:, saccchan, :)))', ...
        squeeze(sum(pupilgrandavg.timelock{sj}(2).lock.trial(:, saccchan, :)))', ...
        squeeze(sum(pupilgrandavg.timelock{sj}(3).lock.trial(:, saccchan, :)))', ...
        squeeze(sum(pupilgrandavg.timelock{sj}(4).lock.trial(:, saccchan, :)))');
    
    % define scalars to do the regression on
    timewin = find(pupilgrandavg.timelock{sj}(4).lock.time < 0 & pupilgrandavg.timelock{sj}(4).lock.time > -0.250 );
    data.blink(find(data.subjnr == sj)) = squeeze(sum(pupilgrandavg.timelock{sj}(4).lock.trial(:, blinkchan, ...
        timewin ), 3));
    data.sacc(find(data.subjnr == sj)) = squeeze(sum(pupilgrandavg.timelock{sj}(4).lock.trial(:, saccchan, ...
        timewin), 3));
    
end

% color scheme
cols = linspecer;

% plot
subplot(4,4,1);
ph = boundedline(1:size(alltimelock.blink, 3), squeeze(nanmean(alltimelock.blink)), ...
    squeeze(nanstd(alltimelock.blink)) / sqrt(length(subjects)), ...
    'cmap', cols);
ylabel({'Blinks'}); axis tight;

% subfunction to put lines and xlabels at the right spots
plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
    pupilgrandavg.timelock{1}(2).lock, [0], ...
    pupilgrandavg.timelock{1}(3).lock, [0 1], ...
    pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
xlabel('Time (ms)');

% plot
subplot(4,4,2);
ph = boundedline(1:size(alltimelock.sacc, 3), squeeze(nanmean(alltimelock.sacc)), ...
    squeeze(nanstd(alltimelock.sacc)) / sqrt(length(subjects)), ...
    'cmap', cols);
ylabel({'Saccades'}); axis tight;
% subfunction to put lines and xlabels at the right spots
plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
    pupilgrandavg.timelock{1}(2).lock, [0], ...
    pupilgrandavg.timelock{1}(3).lock, [0 1], ...
    pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
xlabel('Time (ms)');

print(gcf, '-dpdf', '~/Dropbox/Figures/uncertainty/BlinkSaccadeTimecourse.pdf');


% ==================================================================
% split by uncertainty
% ==================================================================

flds = {'blink', 'sacc'};
% can try this also with all subjects
subjects = 1:27;
figure;
set(gcf, 'DefaultAxesFontSize', 8, 'DefaultAxesFontName', 'Helvetica');

% use nice shades of red and green
cols = linspecer(3); cols = cols(2:3, :);

cnt = 1;
for f = 1:length(flds),
    
    nbins               = 6; % bin in 5 to have comparable plots to the difficulty version?
    data.xval           = abs(data.motionstrength);
    data.rpebin         = nan(size(data.xval)); % preallocate
    
    grandavg.rt.data = nan(length(subjects), 2, nbins);
    
    for sj = subjects,
        
        % loop over error and correct
        cors = [0 1];
        for corr = cors,
            
            % FIT BETAS ON THE FULL MODEL, NOT BINNED
            trls = find(data.subjnr == sj & data.correct == corr);
            
            % dont include RT as a regressor?
            mdl = fitlm(([abs(data.motionstrength(trls))]),  ...
                (data.(flds{f})(trls)));
            
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
                divideintobins(data.xval(trls), data.(flds{f})(trls), nbins);
        end
    end
    
    % PLOT
    subplot(4,4,cnt); cnt = cnt + 1;
    
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
    switch f
        case 1
    ylabel('Blinks');
        case 2
            ylabel('Saccades');
    end
    xlim([0 3.5]); set(gca, 'xtick', 0:1.75:3.5, 'xticklabel', {'hard', 'medium', 'easy' });
    axis square;
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    
    %% make the subplot next to it show the significance of the intercepts
    % and slopes
    
    subplot(4,4,cnt); cnt = cnt + 3;
    
    % slopes
    slopes         = [grandavg.rt.regline(:, 1, 2) grandavg.rt.regline(:, 2, 2)];
    [~, pvalE_interc, ~, stat] = ttest(slopes(:, 1), 0, 'tail', 'both');
    [~, pvalC_interc, ~, stat] = ttest(slopes(:, 2), 0, 'tail', 'both');
    [~, pvalD_interc, ~, stat] = ttest(slopes(:,1), slopes(:,2));
    
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
    sigstar2({[1 2]}, pvalD_interc);
    sigstar2({[1,1], [2,2]}, [pvalE_interc pvalC_interc]);
    ylims = get(gca, 'ylim');  ylim([-max(abs(ylims)) max(abs(ylims))]);
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    %ylim([-0.05 0.1]);
    
end

figpath = '~/Dropbox/Figures/uncertainty';
print(gcf, '-dpdf', sprintf('%s/BlinkSaccadeUncertainty.pdf', figpath));

end


% ==================================================================
% layout, plots lines to indicate event onset
% ==================================================================
function [] = plotLines(refdata, reftp, stimdata, stimtp, respdata, resptp, fbdata, fbtp)

xticks = []; xlabels = {};
for t = 1:length(reftp),
    xticks = [xticks dsearchn(refdata.time', reftp(t))];
    if reftp(t) == 0,
        xlabels = [xlabels 'Stimulus 1'];
    else
        xlabels = [xlabels reftp(t) * 1000];
    end
end

for t = 1:length(stimtp),
    xticks = [xticks length(refdata.time) + dsearchn(stimdata.time', stimtp(t))];
    if stimtp(t) == 0,
        xlabels = [xlabels 'Stimulus 2'];
    else
        xlabels = [xlabels stimtp(t)* 1000];
    end
end

for t = 1:length(resptp),
    xticks = [xticks length(refdata.time) + length(stimdata.time) + ...
        dsearchn(respdata.time', resptp(t))];
    if resptp(t) == 0,
        xlabels = [xlabels 'Response'];
    else
        xlabels = [xlabels resptp(t)* 1000];
    end
end

for t = 1:length(fbtp),
    xticks = [xticks length(refdata.time) + length(stimdata.time) + ...
        length(respdata.time) + dsearchn(fbdata.time', fbtp(t))];
    if fbtp(t) == 0,
        xlabels = [xlabels 'Feedback'];
    else
        xlabels = [xlabels fbtp(t)* 1000];
    end
end

set(gca, 'XTick', xticks, 'XTickLabel', [], ...
    'XTickLabelRotation', -45, 'tickdir', 'out', 'box', 'off');

% add white lines to indicate transitions between intervals
x = length(refdata.time)+.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);
x = length(refdata.time) + length(stimdata.time) +.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);
x = length(refdata.time) + length(stimdata.time) + length(respdata.time) +.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);

% add dotted  black lines to indicate event onset
x = dsearchn(refdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.1);
x = length(refdata.time) + dsearchn(stimdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.1);
x = length(refdata.time) + length(stimdata.time) + dsearchn(respdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k','LineStyle', '-', 'LineWidth', 0.1);
x = length(refdata.time) + + length(stimdata.time) + ...
    length(respdata.time) + dsearchn(fbdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.1);

end
