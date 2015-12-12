% figure 4 overview

close; clear;
lagGroups = 1;
whichmodulator = 'pupil';
fig4a_FruendKernels(lagGroups, whichmodulator);
%fig4b_neuroModDecay(lagGroups);
fig4c_decisionStrategies(lagGroups, whichmodulator);
fig4d_psychFuncShift_Bias(lagGroups, whichmodulator);
fig4e_psychFuncShift_Slope(lagGroups, whichmodulator); 
fig4f_HistoryPupil_Bar(lagGroups, whichmodulator);
fig4g_SjCorrelation(lagGroups, whichmodulator);

suplabel(sprintf('%s, Lags %d %d %d', whichmodulator, lagGroups), 't');
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figure4_%s_lags%d%d%d.pdf', whichmodulator, lagGroups));

