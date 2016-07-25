%% check that this looks the same as RT and pupil for uncertainty encoding
subplot(445); b = Uncertainty_byErrorCorrect('uncertainty');
ylabel('Uncertainty');

% now the repetition behaviour...
nbins = 3;
subplot(446); psychFuncShift_Bias('uncertainty', nbins, []);
title('All trials'); ylim([0.45 0.6]);
subplot(447); psychFuncShift_Bias('uncertainty', nbins, 1);
title('Correct'); ylim([0.45 0.6]);
subplot(448); b =  psychFuncShift_Bias('uncertainty', nbins, 0);
title('Error'); ylim([0.45 0.6]);
