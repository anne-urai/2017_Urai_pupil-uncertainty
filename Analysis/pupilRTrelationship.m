global mypath

% JW found a quadratic relationship between pupil and RT... how about my
% data?
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
nbins = 15; clear dat;

for sj = unique(data.subjnr)',
    [pupmn, rtmn, pupstd, rtstd] = divideintobins(data.response_pupil(sj == data.subjnr), ...
        data.rtNorm(sj == data.subjnr), nbins, 'spearman', @nanmean);
    [pupmd, rtmd, pupiqr, rtiqr] = divideintobins(data.response_pupil(sj == data.subjnr), ...
        data.rtNorm(sj == data.subjnr), nbins, 'spearman', @nanmedian);
    
    dat.pup(sj, :) = pupmn;
    dat.rt(sj, :)  = rtmd;
    dat.rtcoefvar(sj, :) = rtiqr ./ rtmd;
end
    
subplot(331);
h = ploterr(mean(dat.pup), mean(dat.rt), [], ...
    std(dat.rt) ./ sqrt(length(unique(data.subjnr))), 'ko', 'abshhxy', 0);
set(h(1), 'markerfacecolor', 'k', 'markeredgecolor', 'w', 'markersize', 5);
box off; axis square; xlabel('Pupil response (z)'); ylabel('RT (ms)');
offsetAxes;


subplot(332);
h = ploterr(mean(dat.pup), mean(dat.rtcoefvar), [], ...
    std(dat.rtcoefvar) ./ sqrt(length(unique(data.subjnr))), 'ko', 'abshhxy', 0);
set(h(1), 'markerfacecolor', 'k', 'markeredgecolor', 'w', 'markersize', 5);
box off; axis square; xlabel('Pupil response (z)'); ylabel('RT coefvar');
offsetAxes;

print(gcf, '-dpdf', sprintf('%s/Figures/pupilRT.pdf', mypath));


figure;
for sj = 1:27,
    subplot(5,6, sj);
    plot(dat.pup(sj, :), dat.rt(sj, :), '.-');
    axis square; box off; lsline;
end