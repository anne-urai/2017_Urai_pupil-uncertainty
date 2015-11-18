function [] = f3ab_PupilTimecourse()
% show 1. the overall timecourse of the pupil
% 2. pupil timecourse, split by correct and error into bins of stimulus
% difficulty

close all; clear; clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE PUPIL GRAND AVERAGE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/Data/UvA_pupil/GrandAverage/pupilgrandaverage.mat');
subjects = 1:27;

% append all the mean timecourses per condition
for sj = unique(subjects),
    thistable = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data_sj%02d.csv', sj));
    
    pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'EyePupil')==1);
    
    % get all timelock
    alltimelock(sj, 1, :) = cat(2, squeeze(nanmean(pupilgrandavg.timelock{sj}(1).lock.trial(:, pupilchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(2).lock.trial(:, pupilchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, :)))', ...
        squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :)))');
end

% color scheme
cols = linspecer;

% plot
subplot(4,4,1);
ph = boundedline(1:size(alltimelock, 3), squeeze(nanmean(alltimelock)), ...
    squeeze(nanstd(alltimelock)) / sqrt(length(subjects)), ...
    'cmap', cols);
ylabel({'Pupil response' ; '(% signal change)'});

axis tight; set(gca, 'ytick', [-1:2:9], 'ylim', [-1 9]);
% subfunction to put lines and xlabels at the right spots
plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
    pupilgrandavg.timelock{1}(2).lock, [0], ...
    pupilgrandavg.timelock{1}(3).lock, [0 1], ...
    pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
xlabel('Time (ms)');
offsetAxes(gca, 0.12, 0);
print(gcf, '-dpdf', '~/Dropbox/Figures/uncertainty/fig2a_pupil.pdf');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPLIT BY CORR VS ERROR AND DIFFICULTY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tabledat = readtable('~/Data/UvA_pupil/CSV/2ifc_data_allsj.csv');
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens

% get all data
load('~/Data/UvA_pupil/GrandAverage/pupilgrandaverage.mat');

% append all the mean timecourses per condition
for sj = unique(subjects),
    thistable = readtable(sprintf('~/Data/UvA_pupil/CSV/2ifc_data_sj%02d.csv', sj));
    
    cors = [0 1];
    cnt  = 0;
    for c = 1:2,
        
        trls = find(thistable.correct == cors(c));
        motionstrengthquantiles = quantile(abs(thistable.motionstrength(trls)), 2);
        
        for d = 1:3, % divide into 3 bins of absolute motionstrength
            
            % find those trials
            switch d
                case 1
                    trls = find(thistable.correct == cors(c) & abs(thistable.motionstrength) < motionstrengthquantiles(1));
                case 2
                    trls = find(thistable.correct == cors(c) & ...
                        abs(thistable.motionstrength) < motionstrengthquantiles(2) ...
                        &    abs(thistable.motionstrength) > motionstrengthquantiles(1) );
                case 3
                    trls = find(thistable.correct == cors(c) & ...
                        abs(thistable.motionstrength) > motionstrengthquantiles(2) );
            end
            
            pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'EyePupil')==1);
            
            cnt = cnt + 1;
            % get all timelock
            alltimelock(sj, cnt, :) = cat(2, squeeze(nanmean(pupilgrandavg.timelock{sj}(1).lock.trial(trls, pupilchan, :)))', ...
                squeeze(nanmean(pupilgrandavg.timelock{sj}(2).lock.trial(trls, pupilchan, :)))', ...
                squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(trls, pupilchan, :)))', ...
                squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(trls, pupilchan, :)))');
        end
    end
end

% color scheme
cols = cbrewer('div', 'RdYlGn', 10);
cols = cols([1:3 end-2:end], :);

% plot
subplot(4,4,1);
ph = boundedline(1:size(alltimelock, 3), squeeze(nanmean(alltimelock)), ...
    permute(squeeze(nanstd(alltimelock)) / sqrt(length(subjects)), [2 3 1]), ...
    'cmap', cols);
ylabel({'Pupil response' ; '(% signal change)'});
axis tight;
ph2 = plot(1:6, mean(get(gca, 'ylim'))*ones(6, 10), '.w');
lh = legend(ph2); % make t
lh.String = {'\color[rgb]{0.647058823529412,0,0.149019607843137} error hard', ...
    '\color[rgb]{0.843137254901961,0.188235294117647,0.152941176470588} error medium', ...
    '\color[rgb]{0.956862745098039,0.427450980392157,0.262745098039216} error easy', ...
    '\color[rgb]{0.400000000000000,0.741176470588235,0.388235294117647} correct hard', ...
    '\color[rgb]{0.101960784313725,0.596078431372549,0.313725490196078} correct medium', ...
    '\color[rgb]{0,0.407843137254902,0.215686274509804} correct easy'};

lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .15;
set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 6);

axis tight; set(gca, 'ytick', [-1:2:9], 'ylim', [-1 9]);
% subfunction to put lines and xlabels at the right spots
plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
    pupilgrandavg.timelock{1}(2).lock, [0], ...
    pupilgrandavg.timelock{1}(3).lock, [0 1], ...
    pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
xlabel('Time (ms)');
offsetAxes(gca, 0.12, 0);
print(gcf, '-dpdf', '~/Dropbox/Figures/uncertainty/sfig_pupilsplit.pdf');

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

set(gca, 'XTick', xticks, 'XTickLabel', xlabels, ...
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
