%% instead of csv files, write text files that are in the format the Fruend toolbox needs

global mypath;
% cd(sprintf('%s/Code/serial-dependencies/data/', mypath));
cd('~/Dropbox/code/pupilUncertainty/serial-dependencies/');
subjects = 1:27;

for sj = subjects,
    disp(sj);
    % use files with cleaner pupil data
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % remove first session
    data = data(find(data.sessionnr > 1), :);
    
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

    % no modulatrion
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0)];
    
    dlmwrite(sprintf('2ifc_plain_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
    % use this one to add the pupil values
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) zscore(data.decision_pupil)];
    
    dlmwrite(sprintf('2ifc_pupil_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
    % feedback pupil
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) zscore(data.feedback_pupil)];
    
    dlmwrite(sprintf('2ifc_feedbackpupil_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
    
    % feedback pupil with decision pupil regressed out
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) ...
        zscore(projectout(data.feedback_pupil, data.decision_pupil))];
    
    dlmwrite(sprintf('2ifc_fb-decpupil_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
    % double check negative rts
    data.rt(data.rt < 0.01) = 0.01;
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) zscore(log(data.rt))];
    subplot(5,6,sj); histogram(zscore(log(data.rt))); axis tight;
    
    dlmwrite(sprintf('2ifc_rt_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
    % rt projected out of decision pupil
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) ...
        zscore(projectout(data.decision_pupil, data.rt))];
    
    dlmwrite(sprintf('2ifc_pupil-rt_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
     % decision pupil projected out of rt
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) ...
        zscore(projectout(data.rt, data.decision_pupil))];
    
    dlmwrite(sprintf('2ifc_rt-pupil_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
 
    % baseline pupil
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) ...
        zscore(circshift(data.baseline_pupil, -1))];
    
    dlmwrite(sprintf('2ifc_baselinepupil_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
    % evidence strength
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) ...
        abs(data.motionstrength)];
    
    dlmwrite(sprintf('2ifc_evidence_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
    % double modulation: decision pupil and rt
    % log-transform RT, make sure this doesn't include zero!
    newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0) ...
        zscore(data.decision_pupil) zscore(log(data.rt+0.1))];
    dlmwrite(sprintf('2ifc_pupil+rt_sj%02d.txt', sj), ...
        newdat,'delimiter','\t','precision',4);
    
end