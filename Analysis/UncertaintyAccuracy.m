function [] = UncertaintyAccuracy(field)
% plot RT as a function of |motionstrength|, to show that the pattern
% follows model predictions of uncertainty

global mypath;
nbins = 6;

% can try this also with all subjects
subjects = 1:27;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE OVERVIEW OF THE PUPIL UNCERTAINTY CORRELATION FOR ALL THESE DIFFERENT FIELDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grandavg.rt = nan(length(subjects), nbins);
grandavg.acc = nan(length(subjects), nbins);

for sj = subjects,
    
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    [grandavg.rt(sj, :), grandavg.acc(sj, :)] = divideintobins(data.(field), data.correct, nbins);
    
end

% PLOT
h = ploterr(nanmean(grandavg.rt), 100* nanmean(grandavg.acc), nanstd(grandavg.rt) ./ sqrt(length(subjects)), ...
    100 * nanstd(grandavg.acc) ./ sqrt(length(subjects)), 'abshhxy', 0);
set(h(1), 'color', 'k', 'marker', '.', 'markersize', 10);
set(h(2), 'color', 'k');
set(h(3), 'color', 'k');
switch field
    case 'rt'
        xlabel('Reaction time (s)');
        xlim([0.1 1.05])
    case 'decision_pupil'
        xlabel('Pupil response (z)');
        %xlim([0.1 1.05])
end
ylabel('Accuracy (%)');
set(gca, 'ylim', [50 100], 'ytick', 50:25:100);
axis square; box off;

end