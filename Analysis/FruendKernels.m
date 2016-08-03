function [] = FruendKernels(whichmodulator, field)
% history kernels

global mypath;
% determine the subjects based on their plain weights
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator));

if strcmp(whichmodulator, 'plain') && strcmp(field, 'response'),
    % save colors
    colors = cbrewer('div', 'PuOr', 256);
    %colors(1:128, :) = flipud(colors(129:end, :));
    % flip the halves
    colors(128-40:128+40, :) = []; % remove white in the middle
    colors = colors + 0.1; % make a bit paler
    colors(colors > 1) = 1;
    
    for sj = 1:27,
        % find which color this is
        colspace =  linspace(-0.5, 0.5, size(colors, 1));
        colidx   = dsearchn(colspace', dat.(field)(sj, 1));
        mycolmap(sj, :) = colors(colidx, :);
    end
    save(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath), 'mycolmap');
elseif strcmp(field, 'response_pupil'),
    mycolmap = 0.8*ones(27, 3);
else
    % get the colors from all the data combined
    load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));
end

hold on;
for sj = 1:27,
    plot(dat.(field)(sj, :)', 'color', mycolmap(sj, :), 'linewidth', 0.5);
end
scatter(ones(1, 27), dat.(field)(:, 1), 10, mycolmap, 'filled');

if strcmp(whichmodulator, 'plain') && strcmp(field, 'response'),
    
    % show which one is the example
    plot(1, dat.(field)(10, 1), 'ok', 'markersize', 4, 'linewidth', 0.2);
    
    % add the group
    [ax, p1, p2] = plotyy(1:7, nanmean(dat.(field)),  ...
        1:7, mean(abs(dat.(field))));
    
    set(ax(1), 'ycolor', 'k', 'xtick', 1:7, 'ytick', [-0.4 0 0.4], 'ylim', [-.45 .4], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
    set(ax(2), 'ycolor', [0.4 0.4 0.4], 'xtick', 1:7, 'ytick', [0 0.2], 'ylim', [0 .4], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
    
    set(p1, 'color', 'k', 'linewidth', 1);
    set(p2, 'color', [0.4 0.4 0.4], 'linewidth', 1);
    axis(ax(1), 'square'); axis(ax(2), 'square');
    
    ylabel(ax(1),'Choice weight') % label left y-axis
    ylabel(ax(2),'|Choice weight|') % label right y-axis
    
else
    
    plot(1:7, nanmean(dat.(field)), 'k', 'linewidth', 1);
    axis square;
    set(gca, 'ycolor', 'k', 'xtick', 1:7, 'ytick', ...
        [-0.2 0 0.2], 'ylim', [-.25 .25], 'xlim', [0.5 7.5], ...
        'box', 'off', 'xminortick', 'off', 'yminortick', 'off');
    ylabel('Pupil * choice weight');
    
    if strcmp(field, 'response_pupil'),
        % indicate significance for each lag
        clear h;
        for s = 1:size(dat.(field), 2),
            h(s) = ttest(dat.(field)(:, s));
        end
        h(h < 1) = NaN;
        plot(1:7, -0.2*h, 'k.', 'markersize', 5);
    end
end
xlabel('Lags');
end
