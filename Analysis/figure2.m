% replicates figure 2
% Anne Urai, 2015

close; figure;
global mypath;

subplot(441); MotionEnergy_Timecourse;
subplot(442); MotionEnergy_Probability;

% add psychometric functions, chronometric functions and RT distributions
subplot(445); PsychometricFunction;
subplot(446); RTdistributions;
print(gcf, '-dpdf', sprintf('%s/Figures/figure2.pdf', mypath));
