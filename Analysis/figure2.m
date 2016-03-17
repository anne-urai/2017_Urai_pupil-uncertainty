% replicates figure 2
% Anne Urai, 2015

close; figure;
global mypath;

subplot(441); MotionEnergy_Timecourse;
subplot(442); MotionEnergy_Probability;
subplot(445); ReactionTimeUncertainty;

print(gcf, '-dpdf', sprintf('%s/Figures/figure2.pdf', mypath));
