load('~/Data/pupilUncertainty/GrandAverage/pupilgrandaverage.mat');

for sj = 1:27,
    time{sj} = pupilgrandavg.timelock{sj}(1).lock.trialinfo(:, 1);
end

newtime = cat(1, time{:});
difftime = diff(newtime);
difftime(difftime < 0) = [];

difftime = difftime / 100; % divide by sampling rate

% 7 lags
median(difftime)*7;