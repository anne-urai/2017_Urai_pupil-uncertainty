function PsychFuncs_byUncertainty(field)


global mypath;
nbins = 6;

% can try this also with all subjects
subjects = 1:27;
grandavg.ev = nan(length(subjects), 2, nbins);
grandavg.acc = nan(length(subjects), 2, nbins);

for sj = subjects,
    
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % split by low and high RT
    rtMed = quantile(data.(field), 2);
    
    for r = 1:2,
        switch r
            case 1
                trls = find(data.(field) < rtMed(1));
            case 2
                trls = find(data.(field) > rtMed(2));
        end
        [grandavg.ev(sj, r, :), grandavg.acc(sj, r, :)] = ...
            divideintobins(abs(data.motionstrength(trls)), data.correct(trls), nbins);
    end
end

colors(1,:) = [0.8 0.8 0.8];
colors(2,:) = [0.2 0.2 0.2];

% PLOT
hold on;
for r = 1:2,
    h = ploterr(squeeze(nanmean(grandavg.ev(:, r, :))), squeeze(100* nanmean(grandavg.acc(:, r, :))), ...
        squeeze( nanstd(grandavg.ev(:, r, :))) ./ sqrt(length(subjects)), ...
        100 * squeeze(nanstd(grandavg.acc(:, r, :))) ./ sqrt(length(subjects)), 'abshhxy', 0);
    set(h(1), 'color', colors(r,:), 'marker', '.', 'markersize', 10);
    set(h(2), 'color', colors(r,:));
    set(h(3), 'color', colors(r,:));
    handles{r} = h(1);
end

ylabel('Accuracy (%)');
set(gca, 'xlim', [-0.5 5.5], 'ylim', [50 100], 'ytick', 50:25:100);
axis square; box off;
xlabel('Evidence'); xlim([-0.5 5.5]);

switch field
    case 'rt'
        legend([handles{:}], {'Low RT', 'High RT'}, 'location', 'southeast');
        legend boxoff;
    case 'decision_pupil'
        legend([handles{:}], {'Low pupil', 'High pupil'}, 'location', 'southeast');
        legend boxoff;
end
end