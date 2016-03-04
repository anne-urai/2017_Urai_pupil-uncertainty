% visualize designM

sj = 1;
data = readtable(sprintf('~/Data/pupilUncertainty/Data/CSV/2ifc_data_sj%02d.csv', sj));
data = data(1:100, :);
resp = data.resp; resp(resp < 0) = 0;

figure; set(gcf, 'defaultaxesfontsize', 10);
% response
subplot(1,15,1); colormap bone; imagesc(resp); caxis([-1 1]);
set(gca, 'xtick', 1, 'xticklabel', 'responses', 'ytick', [], 'xticklabelrotation', -30);
title('DV');

% basic psychometric function
designM = [ones(100, 1) zscore(data.motionstrength) ...
    circshift(data.resp, 1)  circshift(data.stim, 1) ...
    circshift(zscore(data.decision_pupil), 1)  circshift(zscore(data.decision_pupil), 1).* circshift(data.resp, 1) circshift(zscore(data.decision_pupil), 1).* circshift(data.stim, 1)...
    circshift(zscore(data.rt), 1)  circshift(zscore(data.rt), 1).* circshift(data.resp, 1) circshift(zscore(data.rt), 1).* circshift(data.stim, 1)...
    ];
subplot(1,17,[2 10]); colormap bone; imagesc(designM);
set(gca, 'xtick', 1:size(designM, 2), 'xticklabel', {'intercept', 'slope', 'resp-1', 'stim-1'...
    'pupil-1', 'resp-1*pupil-1', 'stim-1*pupil-1', 'rt-1','resp-1*rt-1', 'stim-1*rt-1'}, ...
    'xticklabelrotation', -30, 'ytick', []);
caxis([-1 1]); title('Design matrix');

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/desigmM.pdf'));
