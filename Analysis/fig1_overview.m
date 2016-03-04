% replicates figure 2
% Anne Urai, 2015

global mypath;

subplot(444); fig1a_PsychometricFunction;
subplot(446); fig1b_MotionEnergy_Timecourse;
subplot(447); fig1b_MotionEnergy_Probability;
subplot(4,4,14); fig1d_ReactionTimeUncertainty;

print(gcf, '-dpdf', sprintf('%s/Figures/figure1.pdf', mypath));
