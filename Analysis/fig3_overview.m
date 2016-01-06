% replicates figure 3
% Anne Urai, 2015

subplot(4,4,1); fig3a_PupilTimecourse;
subplot(4,4,3); fig3b_PupilUncertaintyTimecourse(1);

subplot(4,4,5); fig3c_PupilUncertaintyCorrelation;
subplot(4,4,7); fig3de_Uncertainty_Accuracy;

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty_paper/figure3.pdf'));
