
% reproduces 
global mypath;
close; figure;

% use nice shades of red and green
cols = cbrewer('qual', 'Set1', 8);
cols = cols([1 2], :);

nbins = 3;
subplot(441); psychFuncShift_Bias_Slope('pupil', nbins, 1);
title('Correct', 'color', cols(2,:));
subplot(442); psychFuncShift_Bias_Slope('pupil', nbins, 0);
title('Error', 'color', cols(1,:));

subplot(443); psychFuncShift_Bias_Slope('rt', nbins, 1);
title('Correct', 'color', cols(2,:));
subplot(444); psychFuncShift_Bias_Slope('rt', nbins, 0);
title('Error', 'color', cols(1,:));

print(gcf, '-dpdf', sprintf('%s/Figures/figureS2.pdf', mypath));