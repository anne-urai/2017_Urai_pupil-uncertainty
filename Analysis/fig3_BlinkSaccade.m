function fig3_BlinkSaccade

addpath('~/Documents/fieldtrip');
ft_defaults;
subjects = 1:27;

load('~/Data/pupilUncertainty/GrandAverage/pupilgrandaverage.mat');

% append all the mean timecourses per condition
for sj = unique(subjects),
    blinkchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'Blinks')==1);
    
    % get all timelock
    alltimelock.blink(sj, 1, :) = cat(2, squeeze(nanmean(pupilgrandavg.timelock{sj}(1).lock.trial(:, blinkchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(2).lock.trial(:, blinkchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(:, blinkchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, blinkchan, :)))');
    
    saccchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'Saccades')==1);
    
    % get all timelock
    
    alltimelock.sacc(sj, 1, :) = cat(2, squeeze(nanmean(pupilgrandavg.timelock{sj}(1).lock.trial(:, saccchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(2).lock.trial(:, saccchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(:, saccchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, saccchan, :)))');
    
    % define scalars to do the regression on
    grandavg.blinkcnt(sj, :) = 
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
