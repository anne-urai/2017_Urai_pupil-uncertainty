
close; figure;
global mypath;

%% add psychometric functions, chronometric functions
for s = 2:6,
    subplot(5,5,s-1);
    PsychometricFunction(s);
end

%% run history model for every session
cd(sprintf('%s/Code/serial-dependencies/data', mypath));
subjects = 1:27;

% for supplement, do this per session
for sj = subjects,
    disp(sj);
    
    for session = 2:6,
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        
        if sj == 15,
            % this person did half a session 3 and half a session 6 + full session 7
            % change this so that each chunk has 10 blocks in it
            data.sessionnr(466:937) = 3;
            data.sessionnr(938:1425) = 4;
            data.sessionnr(1426:1921) = 5;
            data.sessionnr(1922:end) = 6;
        end
        
        data = data(find(data.sessionnr == session), :);

        % generate block nrs, NOT identical to session nrs! History effects
        % should not continue beyond a block
        blockchange = find(diff(data.trialnr) < 0);
        blocknrs = zeros(height(data), 1);
        for b = 1:length(blockchange)-1,
            blocknrs(blockchange(b)+1:blockchange(b+1)) = blocknrs(blockchange(b))+1;
        end
        blocknrs(blockchange(end)+1:end) = blocknrs(blockchange(end))+1;
        
        % no modulation, just history
        newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0)];
        
        dlmwrite(sprintf('2ifc_plain_session%d_sj%02d.txt', session, sj), ...
            newdat,'delimiter','\t','precision',4);
    end
    
end

%%

cd(sprintf('%s/Code/serial-dependencies', mypath));
for session = 2:6,
    system([sprintf('for sj in {1..27}; do filename=$(printf "data/2ifc_plain_session%d_sj%%02d.txt" $sj);', session), ...
        sprintf('echo $filename; python2.7 analysis.py -fr -n10 -p "%s/Data/serialmodel" $filename; sleep 5; done', mypath)]);
end

% retrieve those in matlab
cd(sprintf('%s/Code/Analysis', mypath));
for session = 2:6,
    a6_retrieveDataFromPython(sprintf('plain_session%d', session)); % outputs absolute errors
end

print(gcf, '-dpdf', sprintf('%s/Figures/figureS4.pdf', mypath));
