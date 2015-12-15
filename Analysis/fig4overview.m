% figure 4 overview

mods = {'rt', 'pupil', 'feedbackpupil', 'fb-decpupil'};
lags = {[1], [2], [3], [1 2], [1 2 3]};

lags = {[1]};
mods = {'pupil'};

for m = 1:length(mods),
    whichmodulator = mods{m};
    for l = 1:length(lags),
        lagGroups = lags{l};
        
        close; 
        % top row pure history effects
        subplot(441); fig4a_FruendKernels(lagGroups, 'plain', 'individual');
        subplot(442); fig4c_decisionStrategies(lagGroups, 'plain');
        
        subplot(445); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all');
        subplot(446); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'all');
        subplot(447); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'all');
        
        subplot(449); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'repeat');
        subplot(4,4,10); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'repeat');
        subplot(4,4,11); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'repeat');
        
        subplot(4,4,13); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'switch');
        subplot(4,4,14); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'switch');
        
        subplot(4,4,15); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'switch');
        
        suplabel(sprintf('%s, Lags %d %d %d', whichmodulator, lagGroups), 't');
        print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figure4_%s_lags%d%d%d.pdf', whichmodulator, lagGroups));
        
        %% now the supplementary figure
        close all; clc;
        fig4d_psychFuncShift_Bias_byResp(lagGroups, whichmodulator); % split by response
        subplot(4,4,5); fig4g_SjCorrelation(lagGroups, whichmodulator);
        subplot(446); fig4c_decisionStrategies_interaction(lagGroups, whichmodulator);
        subplot(447); fig4S_FruendPlainVsMod(lagGroups, whichmodulator);
        
        suplabel(sprintf('%s, Lags %d %d %d', whichmodulator, lagGroups), 't');
        print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figureS4_%s_lags%d%d%d.pdf', whichmodulator, lagGroups));
    end
end
