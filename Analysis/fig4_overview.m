% figure 4 overview

lagGroups = 1;
mods = {'rt', 'pupil-rt', 'pupil'} %, 'feedbackpupil', 'fb-decpupil'};
%mods = {'evidence'};
nbins = 3;

for m = 1:length(mods),
    whichmodulator = mods{m};
    close;
    figure;
    
    % top row pure history effects
    subplot(441); fig4a_FruendKernels(lagGroups, 'plain', 'individual');
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    subplot(443); fig4b_decisionStrategies(lagGroups, 'plain');
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    % middle row
    subplot(445); fig4c_psychFuncShift_Bias_byResp(lagGroups, whichmodulator, 'all', nbins); % todo
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    subplot(446); fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', nbins, 1);
    hold on; fig4d_psychFuncShift_Bias(lagGroups, whichmodulator, 'all', nbins, 0);
    
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    subplot(447); fig4d_psychFuncShift_Slope(lagGroups, whichmodulator, 'all', [], nbins);
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    
    % bottom row
    subplot(4,4,9); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'all');
    subplot(4,4,10); fig4g_SjCorrelation(lagGroups, whichmodulator);
    
    subplot(4,4,11); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'switch');
    title('Alternators');
    subplot(4,4,12); fig4f_HistoryPupil_Bar(lagGroups, whichmodulator, 'repeat');
    title('Repeaters');
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

    print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figure4_layout2_%s_lags%d%d%d.pdf', whichmodulator, lagGroups));
end
