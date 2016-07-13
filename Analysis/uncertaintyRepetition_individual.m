global mypath;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

% add the code folder
addpath([mypath '/Code/Analysis/']);
cd([mypath '/Code/Analysis/']);
close all;

%% this looks super weird... redo for every individual!
close all; nbins = 4;

for sj = 1:27,
    clf;
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % first, establish the correlation between pupil, RT and uncertainty
    subplot(441); hold on;
    plot(data.uncertainty(data.correct == 1), ...
        data.decision_pupil(data.correct == 1), 'g.');
    plot(data.uncertainty(data.correct == 0), ...
        data.decision_pupil(data.correct == 0), 'r.');
    axis tight; axis square; box off;
    l = lsline;
    l(1).Color = l(1).Color * 0.8;
    l(2).Color = l(2).Color * 0.8;
    
    rho1 = corr(data.uncertainty(data.correct == 1), ...
        data.decision_pupil(data.correct == 1), 'type', 'spearman');
    rho2 = corr(data.uncertainty(data.correct == 0), ...
        data.decision_pupil(data.correct == 0), 'type', 'spearman');
    title(sprintf('%.2f / %.2f', rho1, rho2));
    xlabel('Uncertainty'); ylabel('Pupil');
    
    subplot(442); hold on;
    plot(data.uncertainty(data.correct == 1), ...
        data.rt(data.correct == 1), 'g.');
    plot(data.uncertainty(data.correct == 0), ...
        data.rt(data.correct == 0), 'r.');
    axis tight; axis square; box off;
    l = lsline;
    l(1).Color = l(1).Color * 0.8;
    l(2).Color = l(2).Color * 0.8;
    
    rho1 = corr(data.uncertainty(data.correct == 1), ...
        data.rt(data.correct == 1), 'type', 'spearman');
    rho2 = corr(data.uncertainty(data.correct == 0), ...
        data.rt(data.correct == 0), 'type', 'spearman');
    title(sprintf('%.2f / %.2f', rho1, rho2));
    xlabel('Uncertainty'); ylabel('RT');
    
    subplot(443); hold on;
    plot(data.decision_pupil(data.correct == 1), ...
        data.rt(data.correct == 1), 'g.');
    plot(data.decision_pupil(data.correct == 0), ...
        data.rt(data.correct == 0), 'r.');
    axis tight; axis square; box off;
    l = lsline;
    l(1).Color = l(1).Color * 0.8;
    l(2).Color = l(2).Color * 0.8;
    
    rho1 = corr(data.decision_pupil(data.correct == 1), ...
        data.rt(data.correct == 1), 'type', 'spearman');
    rho2 = corr(data.decision_pupil(data.correct == 0), ...
        data.rt(data.correct == 0), 'type', 'spearman');
    title(sprintf('%.2f / %.2f', rho1, rho2));
    xlabel('Pupil'); ylabel('RT');
    
    % then, show for all those 3
    subplot(445); hold on;
    psychFuncShift_Bias('uncertainty', nbins, 1, sj);
    psychFuncShift_Bias('uncertainty', nbins, 0, sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(1,:) = get(gca, 'ylim');
    
    subplot(446); hold on;
    psychFuncShift_Bias('rt', nbins, 1, sj);
    psychFuncShift_Bias('rt', nbins, 0, sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(2,:) = get(gca, 'ylim');
    
    subplot(447); hold on;
    psychFuncShift_Bias('pupil', nbins, 1, sj);
    psychFuncShift_Bias('pupil', nbins, 0, sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(3,:) = get(gca, 'ylim');
    
    % then, show for all those 3
    subplot(4,4,9); hold on;
    psychFuncShift_Bias('uncertainty', nbins, [], sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(4,:) = get(gca, 'ylim');
    
    subplot(4,4,10); hold on;
    psychFuncShift_Bias('rt', nbins, [], sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(5,:) = get(gca, 'ylim');
    
    subplot(4,4,11); hold on;
    psychFuncShift_Bias('pupil', nbins, [], sj);
    set(gca, 'ytick', 0:0.05:1);
    ylims(6,:) = get(gca, 'ylim');
    
    % make sure the axes are the same
    newylims(1) = min(ylims(:));
    newylims(2) = max(ylims(:));
    for i = [5:7 9:11],
        subplot(4,4,i); ylim(newylims);
    end
    
    print(gcf, '-dpdf', sprintf('%s/Figures/P%02d_uncertaintyRepetition.pdf', mypath, sj));
end


