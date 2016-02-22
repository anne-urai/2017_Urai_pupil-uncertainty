% figure 4 overview

global mypath;

close;
figure;
%s/Data/serialmodel
% top row pure history effects
subplot(441); fig3a_FruendKernels('plain');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

subplot(442); fig3b_decisionStrategies('plain');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

% set correctness to 1 to see the same figures but only based on trials
% after a correct response 
correctness = []; % empty means both correct and error trials will be used
nbins = 3;

subplot(445); fig3c_psychFuncShift_Bias_byResp('pupil', nbins, correctness);
subplot(446); fig3d_psychFuncShift_Bias_Slope('pupil', nbins, correctness);

subplot(4,4,7); fig3c_psychFuncShift_Bias_byResp('rt', nbins, correctness); 
subplot(4,4,8); fig3d_psychFuncShift_Bias_Slope('rt', nbins, correctness);

% bottom row - unique variance of pupil and RT
subplot(4,4,9); fig3hi_HistoryPupil_Bar('pupil', 'all');
subplot(4,4,10); fig3hi_HistoryPupil_Bar('rt', 'all');

subplot(4,4,11); fig3j_SjCorrelation('pupil');
subplot(4,4,12); fig3j_SjCorrelation('rt');

print(gcf, '-dpdf', sprintf('%s/Figures/figure3.pdf', mypath));

%% mean time between responses

load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));

for sj = 1:length(pupilgrandavg.timelock),
   respdiff = diff(pupilgrandavg.timelock{sj}(1).lock.trialinfo(:, 9)) ./ 100;
   trldiff = diff(pupilgrandavg.timelock{sj}(1).lock.trialinfo(:, 12));
   respdiff(trldiff ~= 1) = []; % only use the difference between subsequent trials
   timebetweenResp{sj} = respdiff;
end
timebetweenResp = cat(1, timebetweenResp{:});
median(timebetweenResp); % long-tailed distribution, so mean is biased
