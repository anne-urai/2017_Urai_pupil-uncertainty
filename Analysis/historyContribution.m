
% 10. variance explained by history as a function of stimulus strength
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plainCoh'));
% dat.variance.stimuli    = reshape(dat.variance.stimuli, [27 5 5]);
% dat.variance.explained  = reshape(dat.variance.explained, [27 5 5]);

% first, average across sessions for each person
% http://stackoverflow.com/questions/36965637/matlab-compute-average-value-corresponding-to-each-column-number-in-matrix
dat.newvariance.stimuli     = [0.625 1.25 2.5 5 10 20 30];
dat.newvariance.explained   = nan(27, length(dat.newvariance.stimuli));
for sj = 1:27,
    A = [dat.variance.stimuli(sj, :)' dat.variance.explained(sj, :)'];
    [Aunq,~,Aind] = unique(A(:,1));
    B = [Aunq,accumarray(Aind,A(:,2),[],@mean)] * 100; % percentage of motion coherence and of history contribtuion
    dat.newvariance.explained(sj, ismember(dat.newvariance.stimuli, B(:, 1))) = B(:, 2);
end

subplot(4,4,1); 
boundedline(dat.newvariance.stimuli, ...
    squeeze(nanmean(dat.newvariance.explained)), ...
    squeeze(nanstd(dat.newvariance.explained)) ./ sqrt(sum(~isnan(dat.newvariance.explained))), ...
    'cmap', [0 0 0]);

axis square; box off;
title('TEst');
ylim([-8 100]); xlim([-2 30]);
xlabel('motion coherence'); ylabel('History contribution (%)');

print(gcf, '-dpdf', sprintf('%s/Figures/historyContribution.pdf', mypath));