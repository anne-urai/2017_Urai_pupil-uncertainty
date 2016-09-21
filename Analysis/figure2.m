% replicates figure 2
% Anne Urai, 2015

close all; 
global mypath;

subplot(441); motionEnergy_Timecourse;
subplot(442); motionEnergy_Probability;

% add psychometric functions, chronometric functions and RT distributions
subplot(445); psychometricFunction('all');
subplot(446); b = uncertainty_byErrorCorrect('rt', 6);
print(gcf, '-dpdf', sprintf('%s/Figures/Figure2.pdf', mypath));
