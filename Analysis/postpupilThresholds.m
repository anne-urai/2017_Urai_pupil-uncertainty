
global mypath

mods = {'pupil', 'rt'}; cnt = 1; close all;

load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));
mycolmap = cbrewer('div', 'PuOr', 3);

alternators = find(dat.response(:, 1) < 0);
repeaters = find(dat.response(:, 1) > 0);

for m = 1:length(mods),
    
    nbins = 3; 
    
    % get weibull fit after each bin of pupil repsonse
    grandavg = postPupilBehaviour(mods{m}, nbins, []);
    
    %%
    subplot(4,4,cnt); cnt = cnt + 1;
    plotBetasSwarm([grandavg.weibull(:, 1)  grandavg.weibull(:, end)], [0 0 0; 0 0 0]);
    ylabel('Threshold'); set(gca, 'xtick', 1:2, 'xticklabel', {'lowest', 'highest'});
    ylim([0 5]);
    
    subplot(4,4,cnt);cnt = cnt + 1;
    plotBetasSwarm([grandavg.weibull(alternators, 1)  grandavg.weibull(alternators, end)], [mycolmap(1, :); mycolmap(1, :)]);
    ylabel('Threshold'); set(gca, 'xtick', 1:2, 'xticklabel', {'lowest', 'highest'}); title('Alternators');
    ylim([0 5]);
    
    subplot(4,4,cnt); cnt = cnt + 1;
    plotBetasSwarm([grandavg.weibull(repeaters, 1)  grandavg.weibull(repeaters, end)], [mycolmap(3, :); mycolmap(3, :)]);
    ylabel('Threshold'); set(gca, 'xtick', 1:2, 'xticklabel', {'lowest', 'highest'}); title('Repeaters');
    ylim([0 5]);
    
    subplot(4,4,cnt); cnt = cnt + 1;
    boundedline(1:nbins, mean(grandavg.weibull(alternators, :)), std(grandavg.weibull(alternators, :)) ./ sqrt(length(alternators)), ...
        1:nbins, mean(grandavg.weibull(repeaters, :)), std(grandavg.weibull(repeaters, :)) ./ sqrt(length(repeaters)), 'cmap', mycolmap([1 3], :), 'alpha');
    xlim([0.5 nbins+0.5]);
    ylabel('Threshold');
    
    cnt = cnt + 4;
    
    switch mods{m}
        case 'pupil'
            suplabel('Pupil bins', 'x');
        case 'rt'
            suplabel('RT bins', 'x');
    end
end
print(gcf, '-dpdf', sprintf('%s/Figures/postThresholds.pdf', mypath));
