function scatter3dBinned(data)
global mypath

%% 3d divide into bins
colors = cbrewer('qual', 'Set1', 3);
close all; clc;
nbins = 100;
for c = [0 1],
    thisdat   = data((data.correct == c), :);
    x         = thisdat.evidence;
    y         = thisdat.decision_pupil;
    z         = thisdat.repeat;
    disp(nanmean(z));
    
    % divide into bins of X, then plot Y and Z
    [x, y, z] = divideintobins3(x,y,z, nbins);
    
    % best fitting plane, ignore NaNs
    sf = fit([x', y'], z','poly11');
    
    ax(c+1)   = subplot(2,2,c+1);
    p = plot(sf);
    colormap(cbrewer('seq', 'Greys', 12));
    grid on; % make it look nicer
    set(p,'EdgeColor', 'k','FaceAlpha',0.5)
    hold on;
    
    switch c
        case 0
            cmap = cbrewer('seq', 'Reds', nbins);
        case 1
            cmap = cbrewer('seq', 'Blues', nbins);
    end
    % sort by repetition probability
    [val, idx] = sort(z);
    cmap = cmap(idx, :);
    
    scatter3(x,y,z,40, cmap, 'filled');
    axis tight;
    
    % annotate
    xlabel('evidence');
    ylabel('pupil');
    zlabel('repeat');
    
    switch c
        case 0
            title('Error');
        case 1
            title('Correct');
    end
end

set(findall(gcf,'type','text'),'FontSize',10);
linkprop([ax(1) ax(2)], {'CameraPosition','CameraUpVector'}); %h is the axes handle

end
