subplot(441);
FruendKernels('pupil+rt', 'response_pupil');

%% then separately for the two subgroups
% determine the subjects based on their plain weights
field = 'response_pupil';
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));

for groups = 1:2,
    subplot(4,4,groups+1);

    switch groups
        case 1
            subjects = find(dat.response(:, 1) < 0);
        case 2
            subjects = find(dat.response(:, 1) > 0);
    end
    
    hold on;
    for sj = subjects',
        plot(dat.(field)(sj, :)', 'color', mycolmap(sj, :), 'linewidth', 0.5);
    end
    scatter(ones(1, length(subjects)), dat.(field)(subjects, 1), 10, mycolmap(subjects, :), 'filled');
    
    plot(1:7, nanmean(dat.(field)(subjects, :)), 'k', 'linewidth', 1);
    axis square;
    set(gca, 'ycolor', 'k', 'xtick', 1:7, 'ytick', ...
        [-0.2 0 0.2], 'ylim', [-.25 .25], 'xlim', [0.5 7.5], ...
        'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
    ylabel('Pupil * choice weight');
    
    % indicate significance for each lag
    clear h;
    for s = 1:size(dat.(field), 2),
        h(s) = ttest(dat.(field)(subjects, s));
    end
    h(h < 1) = NaN;
    plot(1:7, -0.2*h, 'k.', 'markersize', 5);
    xlabel('Lags');
    
    switch groups
        case 1
         title('Alternators');
        case 2
         title('Repeaters');
    end
end

print(gcf, '-dpdf', sprintf('%s/Figures/pupilLags.pdf', mypath));
