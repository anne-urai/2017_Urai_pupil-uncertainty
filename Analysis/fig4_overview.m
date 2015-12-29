% figure 4 overview

lags = {[1]};
mods = {'pupil'}; %, 'pupil-rt'};

for m = 1:length(mods),
    whichmodulator = mods{m};
    for l = 1:length(lags),
        lagGroups = lags{l};
        
        close;
        % top row pure history effects
        subplot(441); fig4a_FruendKernels(lagGroups, 'plain', 'individual');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(442); fig4b_decisionStrategies(lagGroups, 'plain');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(443); fig4c_psychFuncShift_Bias_byResp(lagGroups, whichmodulator, 'all'); % todo
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(445); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', 1);
        subplot(445); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(446); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'repeat', 1);
        subplot(4,4,6); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'repeat', 0);
        
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,7); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'switch', 1);
        subplot(4,4,7); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'switch', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        %  suplabel('Pupil response bin on current trial', 'x');
        
        subplot(449); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', 1);
        subplot(449); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,10); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'repeat', 1);
        subplot(4,4,10); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'repeat', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,11); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'switch', 1);
        subplot(4,4,11); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'switch', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        %  suplabel('Pupil response bin on current trial', 'x');
        
        subplot(4,4,13); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'all');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,14); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'repeat');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,15); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'switch');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        suplabel(sprintf('%s, Lags %d %d %d', whichmodulator, lagGroups), 't');
        print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figure4_layout2_%s_lags%d%d%d.pdf', whichmodulator, lagGroups));
        
        if 0,
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
end
