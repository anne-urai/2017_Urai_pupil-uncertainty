
%%
b = UncertaintyAccuracy('decision_pupil');
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));
load(sprintf('%s/Data/GrandAverage/sjcolormap.mat', mypath));
subplot(4,4,1); scatter(dat.response(:, 1), b(:, 2), 5, mycolmap);
ylim([-0.25 0.1]);
[rho, pval] = corr(dat.response(:, 1), b(:, 2));
bf10 = corrbf(rho,27);
axis square;
xlabel('Choice weight'); ylabel('Pupil accuracy weight');
print(gcf, '-dpdf', sprintf('%s/Figures/pupilAccCorrChoiceHistory.pdf', mypath));

%% fit threshold on all data for each SJ
for sj  = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    pBest = fminsearchbnd(@(p) cumWB_LL(p, ...
        abs(data.motionstrength), data.correct), ...
        [1 3 0.1], [0 0 0], [3 6 1]);
    grandavg.b(sj, :) = pBest;
end

load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));
subplot(4,4,1); scatter(dat.response(:, 1), grandavg.b(:, 2), 5, mycolmap);
[rho, pval] = corr(dat.response(:, 1), grandavg.b(:, 2));
bf10 = corrbf(rho,27);
axis square;
fprintf('rho = %.3f, p = %.3f, bf10 = %.3f \n', rho, pval, bf10);
xlabel('Choice weight'); ylabel('Psychometric threshold');
print(gcf, '-dpdf', sprintf('%s/Figures/WbThresholdCorrChoiceHistory.pdf', mypath));

%% sequential biases in stimuli
stimRep.one = nan(27, 7);
stimRep.two = nan(27, 7);

for sj  = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    for s = unique(data.sessionnr)',
        stim1 = abs(data.motionstrength(find(data.sessionnr == s)));
        stim1(diff(data.trialnr(find(data.sessionnr == s))) ~= 1) = NaN;
        stimRep.one(sj, s) = nanmean(abs(diff(stim1)));
        
        stim2 = (data.motionstrength(find(data.sessionnr == s)));
        stim2(diff(data.trialnr(find(data.sessionnr == s))) ~= 1) = NaN;
        stimRep.two(sj, s) = nanmean(diff(stim2));
    end
end

subplot(221);
histogram(stimRep.one(:));
box off;
xlabel('Sequential motion energy');

subplot(222);
histogram(stimRep.two(:));
