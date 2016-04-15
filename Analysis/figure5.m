
% reproduces 
global mypath;
close; figure;

% set correctness to 1 to see the same figures but only based on trials
% after a correct response 
correctness = []; % empty; both correct and error trials will be used
% change to 1 or 0 to see the data and get stats for repetition bias after
% correct or error trials, respectively.
nbins = 3;

subplot(441); psychFuncShift_Bias_byResp('pupil', nbins, correctness);
subplot(442); psychFuncShift_Bias_Slope('pupil', nbins, correctness);

subplot(443); psychFuncShift_Bias_byResp('rt', nbins, correctness); 
subplot(444); psychFuncShift_Bias_Slope('rt', nbins, correctness);

print(gcf, '-dpdf', sprintf('%s/Figures/figure5.pdf', mypath));

%%  redo for evidence strength
nbins = 3;

subplot(441); psychFuncShift_Bias_Slope('evidence', nbins, 1);
title('Correct trials');
subplot(442); psychFuncShift_Bias_Slope('evidence', nbins, 0);
title('Error trials');

subplot(443); psychFuncShift_Bias_Slope('pupil', nbins, 1);
title('Correct trials');
subplot(444); psychFuncShift_Bias_Slope('pupil', nbins, 0);
title('Error trials');

