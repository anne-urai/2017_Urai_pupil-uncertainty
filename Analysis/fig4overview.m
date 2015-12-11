% figure 4 overview

close; clear;
lagGroups = 1:2;
fig4a_FruendKernels(lagGroups);
fig4b_neuroModDecay(lagGroups);
fig4c_decisionStrategies(lagGroups);
fig4d_psychFuncShift_Bias(lagGroups);
fig4e_psychFuncShift_Slope(lagGroups); % to do: adapt
fig4f_HistoryPupil_Bar(lagGroups);

suplabel(sprintf('Lags %d %d %d', lagGroups), 't');
print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figure4_lags%d%d%d.pdf', lagGroups));
