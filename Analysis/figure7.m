% figure 4 overview
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'pupil+rt'));
subplot(4,4,1); plotBetasSwarm([dat.response_pupil(:, 1) ...
    dat.stimulus_pupil(:, 1)  dat.response_rt(:, 1) dat.stimulus_rt(:, 1)], ...
    [0 0 0; 0 0 0; 0 0 0; 0 0 0]);
set(gca, 'xtick', 1:4, 'xticklabel', ...
    {'Pupil * choice', 'Pupil * stimulus', 'RT * choice', 'RT * stimulus'}, ...
    'xticklabelrotation', -30);

subplot(4,4,3); SjCorrelation('pupil', 'response');
subplot(4,4,4); SjCorrelation('rt', 'response');

%% show median split for correlation stuff
subplot(4,6, 13); MedianSplit('pupil', 'response'); 
ylim([-0.3 0.15]); ylabel('Pupil * choice');
subplot(4,6,14); MedianSplit('rt', 'response'); 
ylim([-0.3 0.15]); set(gca, 'yaxislocation', 'right');
ylabel('RT * choice');

% group split, also correct and error
subplot(4,6, 15); MedianSplit('pupil', 'correct');
ylim([-0.35 0.2]);
ylabel('Pupil * correct');
subplot(4,6,16); MedianSplit('pupil', 'incorrect');
set(gca, 'yaxislocation', 'right');
ylabel('Pupil * error');ylim([-0.35 0.2]);

% group split, also correct and error
subplot(4,6, 17); MedianSplit('rt', 'correct');
ylim([-0.35 0.2]);
ylabel('RT * correct');
subplot(4,6,18); MedianSplit('rt', 'incorrect');
set(gca, 'yaxislocation', 'right');
ylabel('RT * error');ylim([-0.35 0.2]);

print(gcf, '-dpdf', sprintf('%s/Figures/figure7.pdf', mypath));

%% 
subplot(4,6,19);  SjCorrelation('pupil', 'correct');
subplot(4,6,20);  SjCorrelation('pupil', 'incorrect');




%% mean time between responses

load(sprintf('%s/Data/GrandAverage/pupilgrandaverage.mat', mypath));
timebewteenResp = {};
for sj = 1:length(pupilgrandavg.timelock),
   respdiff = diff(pupilgrandavg.timelock{sj}(1).lock.trialinfo(:, 9)) ./ 100;
   trldiff = diff(pupilgrandavg.timelock{sj}(1).lock.trialinfo(:, 12));
   respdiff(trldiff ~= 1) = []; % only use the difference between subsequent trials
   timebetweenResp{sj} = respdiff;
end
timebetweenResp = cat(1, timebetweenResp{:});
median(timebetweenResp); % long-tailed distribution, so mean is biased


%% check that bias or RT does not change following uncertainty
biasRT;
% these stats are reported in the text, data not shown.

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

% correlate with fruend stimulus and response weights
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));

[rho, pval] = corr(meantransprob', dat.response(:, 1))
bf10 = corrbf(rho,27);
disp(1/bf10);

[rho, pval] = corr(meantransprob', dat.stimulus(:, 1))
bf10 = corrbf(rho,27);
disp(1/bf10)

%%
close;
for sj = 1:27,
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
[newx, newy, stdx, stdy] = divideintobins(abs(data.motionstrength), zscore(data.decision_pupil), 20);
subplot(5,6,sj); plot(newx, newy, 'o'); axis tight; box off;
lsline; 
end

