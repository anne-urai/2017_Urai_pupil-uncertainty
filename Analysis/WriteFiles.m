%% WriteFiles
cd ~/code/serial-dependencies/data

figure; clear;
subjects = 1:27;
for sj = subjects,
    disp(sj);
    % use files with cleaner pupil data
    data = readtable(sprintf('~/Data/HD1/UvA_pupil/CSV/2ifc_data2_sj%02d.csv', sj));
    
    % generate block nrs, NOT identical to session nrs! History effects
    % should not continue beyond a block
    blockchange = find(diff(data.trialnr) < 0);
    blocknrs = zeros(height(data), 1);
    for b = 1:length(blockchange)-1,
        blocknrs(blockchange(b)+1:blockchange(b+1)) = blocknrs(blockchange(b))+1;
    end
    blocknrs(blockchange(end)+1:end) = blocknrs(blockchange(end))+1;
    
    % check how this looks
    % plot(blocknrs, table2array(data(:, 10)));
    % waitforbuttonpress;
    
    % 16.11, use abs(motionstrength) because the toolbox will multiply with
    % stim identity again
    % newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.stim > 0) (data.resp > 0)];
    
    % use this one to add the pupil values
    newdat = [ blocknrs data.sessionnr abs(data.coherence) (data.stim > 0) (data.resp > 0) zscore(log(data.rt))];
    subplot(5,6,sj); histogram(zscore(log(data.rt)), 'binwidth', 0.001);

    
    dlmwrite(sprintf('2ifc_rt_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
end