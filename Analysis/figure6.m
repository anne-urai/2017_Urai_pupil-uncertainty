% figure 4 overview


% bottom row - unique variance of pupil and RT
subplot(4,4,9); fig3h_HistoryPupil_Bar('pupil', 'all');
subplot(4,4,10); fig3h_HistoryPupil_Bar('rt', 'all');

subplot(4,4,11); fig3i_SjCorrelation('pupil');
subplot(4,4,12); fig3i_SjCorrelation('rt');

% show median split for correlation stuff
subplot(4,4,14); fig3j_MedianSplit('pupil'); 
subplot(4,4,15); fig3j_MedianSplit('rt'); 


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


%% check that bias or RT does not change following uncertainty
fig3z_biasRT;

%% compute the actual stimulus transition probabilites

for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    trldiff = diff(data.trialnr);
    blockdiff = diff(data.blocknr);
    transitions = (trldiff ~= 1 | blockdiff ~= 0);
    transitions = [0; transitions];
    data.transitions = transitions;
    
    % loop over the data
    transitions = [1; find(data.transitions == 1); length(data.transitions)];
    transitionprobability{sj} = nan(1, length(transitions)-1);
    
    for b = 1:length(transitions)-1,
        thisdat = data(transitions(b) : transitions(b+1)-1, :);
        if height(thisdat) == 50,
            transitionprobability{sj}(b) = mean((abs(diff(thisdat.stim)) > 0));
        end
    end
    meantransprob(sj) = nanmean(transitionprobability{sj});
end

% correlate with fruend stimulus weights
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));

plot(meantransprob, dat.response(:, 1), 'o');
[rho, pval] = corr(meantransprob', dat.response(:, 1))
bf10 = corrbf(rho,27)

%%
close;
for sj = 1:27,
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
[newx, newy, stdx, stdy] = divideintobins(abs(data.motionstrength), zscore(data.decision_pupil), 20);
subplot(5,6,sj); plot(newx, newy, 'o'); axis tight; box off;
lsline; 
end

