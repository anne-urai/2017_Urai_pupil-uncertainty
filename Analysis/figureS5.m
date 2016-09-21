
clearvars -except mypath;
global mypath; close all; clc;
mods = {'pupil', 'rt'};
nbins = 3;
for m = 1:2,
    subplot(4,4,m);
    grandavg{m} = postPupilBehaviour(mods{m}, 3, []);
    plot([1 nbins], [0.5 0.5], 'k:', 'linewidth', 0.5); hold on;
    
    for r = 1:2,
        x = [1:3] + (r-1)*0.1;
        y = squeeze(grandavg{m}.logisticHistory(:, r, :));
        colors = cbrewer('qual', 'Set2', 3);
        switch r
            case 1
                thismarker = '.';
                thismarkersize = 14;
                thiscolor = colors(1,:);
            case 2
                thismarker = '.';
                thismarkersize = 14;
                thiscolor = colors(3,:);
        end
        
        h = ploterr(x, squeeze(nanmean(y)), [], ...
            squeeze(nanstd(y)) ./sqrt(length(y)), 'k-',  'abshhxy', 0);
        set(h(1), 'color', thiscolor, 'markersize', thismarkersize, 'marker',thismarker);
        set(h(2), 'color', thiscolor); % line color
    end
    
    switch mods{m}
        case 'pupil'
            xticklabs       = repmat({' '}, 1, nbins);
            xticklabs{1}    = 'low';
            xticklabs{end}  = 'high';
            if nbins == 3, xticklabs{2} = 'med'; end
        case 'rt'
            xticklabs       = repmat({' '}, 1, nbins);
            xticklabs{1}    = 'fast';
            xticklabs{end}  = 'slow';
            if nbins == 3, xticklabs{2} = 'med'; end
    end
    
    set(gca, 'xlim', [0.5 nbins+0.5], 'xtick', 1:nbins,  'xticklabel', xticklabs, ...
        'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5, 'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
    axis square; xlim([0.5 nbins+0.5]);
    ylabel('P(choice = 1');
    
    % do statistics, repeated measures anova across bins
    switch mods{m}
        case 'pupil'
            xlabel('Previous trial pupil');
            ylim([0.45 0.55]); set(gca, 'ytick', 0.45:0.05:0.55);
        case 'rt'
            xlabel('Previous trial RT');
            ylim([0.43 0.58]); set(gca, 'ytick', 0.43:0.05:0.58);

        otherwise
            xlabel(sprintf('Previous trial %s', mods{m}));
    end
end

print(gcf, '-dpdf', sprintf('%s/Figures/FigureS5.pdf', mypath));
