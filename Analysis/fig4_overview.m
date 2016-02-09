% figure 4 overview

global mypath;

close;
figure;
%s/Data/serialmodel
% top row pure history effects
subplot(441); fig4a_FruendKernels('plain');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

subplot(442); fig4b_decisionStrategies('plain');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

% middle row
nbins = 3;
subplot(445); fig4c_psychFuncShift_Bias_byResp('pupil', nbins);
subplot(446); fig4d_psychFuncShift_Bias_Slope('pupil', nbins);

subplot(4,4,7); fig4c_psychFuncShift_Bias_byResp('rt', nbins); 
subplot(4,4,8); fig4d_psychFuncShift_Bias_Slope('rt', nbins);

% bottom row - unique variance of pupil and RT
subplot(4,4,9); fig4hi_HistoryPupil_Bar('pupil', 'all');
subplot(4,4,10); fig4hi_HistoryPupil_Bar('rt', 'all');

subplot(4,4,11); fig4j_SjCorrelation('pupil');
subplot(4,4,12); fig4j_SjCorrelation('rt');

print(gcf, '-dpdf', sprintf('%s/Figures/figure4.pdf', mypath));


