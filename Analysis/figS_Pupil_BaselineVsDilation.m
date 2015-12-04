%function [] = s3c_Pupil_BaselineVsDilation()
% from the individual pupil data

subjects = 1:27;
for whichpup = 1,
    
    close all;
    figure; set(gcf, 'DefaultAxesFontSize', 5, 'color', 'w');
    clc;
    
    switch whichpup
        case 1
            xlab = 'baseline_pupil';
            ylab = 'decision_pupil';
            
        case 2
            xlab = 'decision_pupil';
            ylab = 'fbbaseline_pupil';
            
        case 3
            xlab = 'fbbaseline_pupil';
            ylab = 'feedback_pupil';
            
        case 4
            xlab = 'decision_pupil';
            ylab = 'feedback_pupil';
    end
    
    for sj = subjects,
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PREPARE
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
        
        subplot(7,7,find(sj==subjects));
        plot(data.(xlab), data.(ylab), '.');
        axis tight; l = lsline; set(l, 'color', 'k');
        axis square;
        
        box off;
        offsetAxes;
        
        [rho, pval] = corr(data.(xlab), data.(ylab), 'type', 'Spearman');
        title(sprintf('rho %.2f', rho));
        
        allrho(find(sj==subjects)) = rho;
        pv(find(sj==subjects)) = pval;
        
    end
    
    %suplabel(regexprep(xlab, '_', ' '), 'x'); suplabel(regexprep(ylab, '_', ' '), 'y');
    suplabel('Baseline pupil (% signal change)', 'x');
    suplabel('Decision pupil (% signal change from baseline', 'y');
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, '-dpdf', '-painters', sprintf('~/Data/pupilUncertainty/Figures/pupilBaselineVDilation_%d.pdf', whichpup));
    
    % see if the correlation is different from zero
    disp(whichpup)
    mean(allrho)
    %  [~, pval] = permtest(allrho)
    
    %end
end
