% figure 4 overview

lags = {[1]};
mods = {'pupil'}; %, 'pupil-rt'};
%mods = {'feedbackpupil', 'fb-decpupil'};

for m = 1:length(mods),
    whichmodulator = mods{m};
    for l = 1:length(lags),
        lagGroups = lags{l};
        
        close;
        figure;
        
        % top row pure history effects
        subplot(441); fig4a_FruendKernels(lagGroups, 'plain', 'individual');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(442); fig4b_decisionStrategies(lagGroups, 'plain');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        % subplot(4,4,4); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'all');
        % set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        % middle row
        subplot(445); fig4c_psychFuncShift_Bias_byResp(lagGroups, whichmodulator, 'all'); % todo
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(446); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', []);
        % subplot(446); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(447); fig4_psychFuncShift_RT(lagGroups, whichmodulator, 'all', []);
        % subplot(446); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(448); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', []);
        %  subplot(447); fig4e_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        % bottom row
        subplot(449); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'switch');
        % subplot(449); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'switch', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,4,10); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'repeat');
        % subplot(4,4,10); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'repeat', 0);
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        % bottom row
        % whichmodulator =  'pupil_nostim';
        subplot(4,12,32); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'switch');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        subplot(4,12,33); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'repeat');
        set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
        
        %  subplot(4,4,12);
        %  fig4g_SjCorrelation(lagGroups, whichmodulator)
        
        switch whichmodulator
            case 'baseline'
                % extra, baseline pupil stuff
                subplot(4,4,13);
                fig4g_neuroModDecay(lagGroups, whichmodulator, 'repeat', 1);
                fig4g_neuroModDecay(lagGroups, whichmodulator, 'repeat', 0);
                
                subplot(4,4,14);
                fig4g_neuroModDecay(lagGroups, whichmodulator, 'switch', 1);
                fig4g_neuroModDecay(lagGroups, whichmodulator, 'switch', 0);
        end
        print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figure4_layout2_%s_lags%d%d%d.pdf', whichmodulator, lagGroups));
        
    end
end
