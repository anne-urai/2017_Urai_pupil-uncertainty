% replicates figure 2
% Anne Urai, 2015

subplot(444); fig2a_PsychometricFunction;

subplot(446); fig2b_MotionEnergy_Timecourse;
subplot(447); fig2b_MotionEnergy_Probability;

subplot(4,4,14); fig2d_ReactionTimeUncertainty;

print(gcf, '-dpdf', sprintf('%s/Figures/figure2.pdf', mypath));
