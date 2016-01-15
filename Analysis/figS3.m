% figure 4 overview

lags = {[1]};
mods = {'feedbackpupil', 'fb-decpupil'};

close;
figure;

for m = 1:length(mods),
    whichmodulator = mods{m};
    for l = 1:length(lags),
        lagGroups = lags{l};
        
        % top row pure history effects
        subplot(4,4,(m-1)*8 + 1); fig4c_psychFuncShift_Bias_byResp(lagGroups, whichmodulator, 'all'); % todo
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,(m-1)*8 + 2); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', 1);
        subplot(4,4,(m-1)*8 + 2); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,(m-1)*8 + 3); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', 1);
        subplot(4,4,(m-1)*8 + 3); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,(m-1)*8 + 4); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'all');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        % middle row
         subplot(4,4,(m-1)*8 + 5); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'switch', 1);
        subplot(4,4,(m-1)*8 + 5); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'switch', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
      
        subplot(4,4,(m-1)*8 + 6); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'repeat', 1);
        subplot(4,4,(m-1)*8 + 6); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'repeat', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        % bottom row
        subplot(4,8, 13 + (m-1)*16); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'switch');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        subplot(4,8, 14 + (m-1)*16); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'repeat');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,(m-1)*8 + 8);
        fig4g_SjCorrelation(lagGroups, whichmodulator)
        
    end
    
end
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figureS3.pdf'));
