
% reproduces
global mypath;
close; figure;

%% correct and error weights (not the pupil modulation)
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'))

subplot(441);
hold on;
for sj = 1:27,
    h = ploterr(dat.correct(sj, 1), dat.incorrect(sj, 1), ...
        [], ...
        [], '.', 'abshhxy', 0);
    set(h(1), 'color', mycolmap(sj, :), 'markerfacecolor', mycolmap(sj, :), 'markersize', 10);
end

% layout
xlim([-0.5 0.5]); ylim([-1 1]);
xlabel('Reward weight'); ylabel('Punishment weight');
box on; axis square;

%% also show the bar graphs of the error and correct fruend-derived model
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));

colors = cbrewer('qual', 'Set1', 8);
subplot(4,4,2); plotBetasSwarm([dat.correct_pupil(:, 1) dat.incorrect_pupil(:, 1)], colors([2 1], :)); 
set(gca, 'xtick', 1:2, 'xticklabel', {'pupil*correct', 'pupil*error'}, 'xticklabelrotation', -30); %ylim([-0.35 0.3]);
axis square;

colors = cbrewer('qual', 'Set1', 8);
subplot(4,4,3); plotBetasSwarm([dat.correct_rt(:, 1) dat.incorrect_rt(:, 1)], colors([2 1], :)); 
set(gca, 'xtick', 1:2, 'xticklabel', {'rt*correct', 'rt*error'}, 'xticklabelrotation', -30); %ylim([-0.35 0.3]);
axis square;

%% use nice shades of red and green
cols = cbrewer('qual', 'Set1', 8);
cols = cols([1 2], :);

nbins = 3;
subplot(445); psychFuncShift_Bias('pupil', nbins, 1);
title('Correct', 'color', cols(2,:));
subplot(446); psychFuncShift_Bias('pupil', nbins, 0);
title('Error', 'color', cols(1,:));

subplot(447); psychFuncShift_Bias('rt', nbins, 1);
title('Correct', 'color', cols(2,:));
subplot(448); psychFuncShift_Bias('rt', nbins, 0);
title('Error', 'color', cols(1,:));



%%
print(gcf, '-dpdf', sprintf('%s/Figures/figureS2.pdf', mypath));

%% simulate error versus correct and choice versus stimulus weights
% 
% subjects = 1:27;
% for sj = unique(subjects),
%     
%     data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
%     data = data((data.sessionnr > 1), :); % get rid of residual learning effects
%     
%     % get the variables we need
%     resp = data.resp; resp(resp == -1) = 0; % predict response identity
%     motionstrength = zscore(data.motionstrength);
%     prevResp       = circshift(data.resp, 1);
%     prevStim       = circshift(data.stim, 1);
%     prevReward     = circshift(data.correct, 1) .* circshift(data.resp, 1); % Busse: 
%     prevPunishment = circshift(~data.correct, 1) .* circshift(data.resp, 1);
%     
%     % don't use trials that are at the beginning of each block
%     trlDif = [0; diff(data.trialnr)];
%     removeTrls = false(size(trlDif));
%     removeTrls(trlDif < 1) = true;
%     removeTrls(trlDif > 1) = true;
%     removeTrls(find(trlDif > 1) + 1) = true;
%     
%     % =============================== %
%     % 1. Fruend-style
%     % =============================== %
%     
%     % make design matrix - intercept will be added automatically
%     designM = [motionstrength prevResp prevStim];
%     % in these trials, the history wont be able to predict the response
%     designM(removeTrls==1, 2:end) = 0;
%     
%     % fit
%     mdl = fitglm(designM, resp, ...
%         'distr', 'binomial', 'link', 'logit');
%     grandavg.choice(sj) = mdl.Coefficients.Estimate(3);
%     grandavg.stim(sj) = mdl.Coefficients.Estimate(4);
%     
%     % =============================== %
%     % 2. Busse-style
%     % =============================== %
%     
%     designM = [motionstrength prevReward prevPunishment];
%     % in these trials, the history wont be able to predict the response
%     designM(removeTrls==1, 2:end) = 0;
%     
%     % fit
%     mdl = fitglm(designM, resp, ...
%         'distr', 'binomial', 'link', 'logit');
%     grandavg.reward(sj) = mdl.Coefficients.Estimate(3);
%     grandavg.punishment(sj) = mdl.Coefficients.Estimate(4);
%  
% end
% 
% % =============================== %
% % compute Fruend-style correct and error weights
% % =============================== %
% 
% grandavg.correct = grandavg.choice + grandavg.stim;
% grandavg.error   = grandavg.choice - grandavg.stim;
% 
% % plot
% subplot(441); scatter(grandavg.correct, grandavg.reward); xlabel('Correct'); ylabel('Reward');
% subplot(442); scatter(grandavg.error, grandavg.punishment); xlabel('Error'); ylabel('Punishment');
% 
% %% show the same plot as in Figure 6 but with error and correct
% 
% subplot(4,4,1); SjCorrelation('pupil', 'correct');
% subplot(4,4,2); SjCorrelation('rt', 'correct');
% 
% subplot(4,4,3); SjCorrelation('pupil', 'incorrect');
% subplot(4,4,4); SjCorrelation('rt', 'incorrect');