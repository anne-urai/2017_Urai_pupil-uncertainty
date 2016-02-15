function [] = fig3a_PupilTimecourse(plotTotalAverage)
% show 1. the overall timecourse of the pupil
% 2. pupil timecourse, split by correct and error into bins of stimulus
% difficulty

global mypath;

if plotTotalAverage,
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT THE PUPIL GRAND AVERAGE
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));
    subjects = 1:27;
    
    % append all the mean timecourses per condition
    for sj = unique(subjects),
        pupilchan       = find(strcmp(pupilgrandavg.timelock{sj}(1).lock.label, 'EyePupil')==1);
        
        % get all timelock
        alltimelock(sj, 1, :) = cat(2, squeeze(nanmean(pupilgrandavg.timelock{sj}(1).lock.trial(:, pupilchan, :)))', ...
            squeeze(nanmean(pupilgrandavg.timelock{sj}(2).lock.trial(:, pupilchan, :)))', ...
            squeeze(nanmean(pupilgrandavg.timelock{sj}(3).lock.trial(:, pupilchan, :)))', ...
            squeeze(nanmean(pupilgrandavg.timelock{sj}(4).lock.trial(:, pupilchan, :)))');
    end
    
    % plot
    subplot(4,4,1);
    ph = boundedline(1:size(alltimelock, 3), squeeze(nanmean(alltimelock)), ...
        squeeze(nanstd(alltimelock)) / sqrt(length(subjects)), ...
        'cmap', [0 0 0]);
    ylabel({'Pupil response (z)'});
    
    axis tight; set(gca, 'ytick', [0:0.5:1], 'ylim', [-0.2 1]);
    % subfunction to put lines and xlabels at the right spots
    plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
        pupilgrandavg.timelock{1}(2).lock, [0], ...
        pupilgrandavg.timelock{1}(3).lock, [0 1], ...
        pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
    %xlabel('Time (ms)');
    % offsetAxes(gca, 0.12, 0);
    %print(gcf, '-dpdf', '~/Dropbox/Figures/uncertainty/overallPupilTimecourse.pdf');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPLIT BY CORR VS ERROR AND DIFFICULTY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except mypath plotTotalAverage; subjects = 1:27;
warning('error', 'stats:LinearModel:RankDefDesignMat'); % stop if this happens

% get all data
load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));

% append all the mean timecourses per condition
for sj = unique(subjects),
    thistable = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
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
cols = cols([3 2 1 end-2:end], :);

if plotTotalAverage,
    subplot(4,4,5);
end

% plot
ph = boundedline(1:size(alltimelock, 3), squeeze(nanmean(alltimelock)), ...
    permute(squeeze(nanstd(alltimelock)) / sqrt(length(subjects)), [2 3 1]), ...
    'cmap', cols);
ylabel('Pupil response (z)');
axis tight;
ph2 = plot(1:6, mean(get(gca, 'ylim'))*ones(6, 10), '.w');
lh = legend(ph2); % make t
legnames = {'error weak', 'error medium', 'error strong', 'correct weak', 'correct medium', 'correct strong'};
for i = 1:size(cols, 1),
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', cols(i, 1), cols(i, 2), cols(i, 3), legnames{i})];
end
lh.String = str;

lpos = get(lh, 'Position'); lpos(1) = lpos(1) + .15;
set(lh, 'Position', lpos, 'box', 'off', 'FontSize', 6);

axis tight; set(gca, 'ytick', [0:0.5:1], 'ylim', [-0.2 1]);
% subfunction to put lines and xlabels at the right spots
plotLines(pupilgrandavg.timelock{1}(1).lock, [0], ...
    pupilgrandavg.timelock{1}(2).lock, [0], ...
    pupilgrandavg.timelock{1}(3).lock, [0 1], ...
    pupilgrandavg.timelock{1}(4).lock, [0 1 2]);
% xlabel('Time (ms)');

end

