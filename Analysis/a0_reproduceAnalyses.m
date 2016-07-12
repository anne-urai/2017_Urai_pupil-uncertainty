% Urai AE, Braun A, Donner TH. Decision uncertainty drives pupil-linked
% arousal systems and alters subsequent decision-making.

% first, make sure this path matches the place where the data are stored
% and you unzipped everything
global mypath;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

% add the code folder
addpath([mypath '/Code/Analysis/']);
cd([mypath '/Code/Analysis/']);

% optional: get my Tools folder here https://github.com/anne-urai/Tools and
% add this to your path as well - or manually get those functions that
% you're missing (mainly for plotting and some stats).

% if you're not interested in preprocessing and want to quickly reproduce
% the figures, I'd recommend getting just the CSV, serialmodel and GrandAverage 
% folders that will allow you to skip the first 4 steps described here.

%% process pupil data
% note: you can skip this and download the 'P??_alleye.mat' files as well.

% download and install fieldtrip - see http://www.fieldtriptoolbox.org
addpath('~/Documents/fieldtrip/');
ft_defaults;

for sj = 1:27, a1_PupilAnalysis(sj); end

%% motion energy filtering
% note: since this very is memory and CPU heavy, the outputs of this script
% have been provided in motionEnergy.zip - consider using those to continue
% if you don't have access to a computing cluster.
for sj = 1:27, a2_MotionEnergy(sj); end

%% combine all the data that was processed into one grand average file
% note: the resulting file is also included in the GrandAverage folder
a3_writeData2GA;

%% also write one CSV file per sj for quick behavioural data analysis
% note: these files are also included in CSV.zip, so you could skip
% everything up to here if you're not interested in visualizing the pupil
% timecourses
a4_writeData2CSV;

%% reproduce figure 1, 2, 3
close all;
figure1;

close all;
figure2;

close all;
figure3; 

%% run the python model with modulatory term
% write away text files that have the format Python needs
a5_writeFiles4pythonToolbox;

% then, call the terminal from Matlab (not sure if this would work on Windows)
mods = {'plain', 'pupil+rt', 'fbpupil', 'fb+decpupil'};
if ~exist(sprintf('%s/Data/serialmodel', mypath), 'dir'), mkdir(sprintf('%s/Data/serialmodel', mypath)); end

% this is easiest to run from the terminal. With Python 2.7 installed, go
% to the folder mypath/Code/serial-dependencies
cd(sprintf('%s/Code/serial-dependencies', mypath));

for m = 1:length(mods),
    system([sprintf('for sj in {1..27}; do filename=$(printf "data/2ifc_%s_sj%%02d.txt" $sj);', mods{m}), ...
        sprintf('echo $filename; python2.7 analysis.py -fr -n1000 -p "%s/Data/serialmodel/" $filename; sleep 5; done', mypath)]);
end

% important: if you want the quick and dirty version without accurate
% errorbars, change -n1000 to n-10. Otherwise, the code will run bootstraps
% which can take quite a while.

%% get the output back into something Matlab can work with
cd(sprintf('%s/Code/Analysis', mypath));
for m = 1:length(mods),
    a6_retrieveDataFromPython(mods{m}); % outputs absolute errors
end

%% reproduce figures 4,5,6 based on the model fits
close all;
figure4;

close all;
figure5;

close all; 
figure6;

%% also ccreate the supplementary figures
figureS1_MotionEnergy_Filters;
figureS2_RT;
figureS3_feedbackpupil;
figureS4_History_CorrectError;
figureS5_pupilResponseLagged;
figureS5_scatterIndividual;
figureS6_performanceOverSessions;

%% analyze some final stuff that's not in the figures

% 1. median time between responses
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

% 2.  Speed-accuracy trade-off 
biasRT;

% 3.  compute autocorrelation in evidence strength
stimRep.rho = nan(27, 7);
stimRep.pval = nan(27, 7);
for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    for session = unique(data.sessionnr)',
        thisdat = data(find(session==data.sessionnr), :);
        stim1 = abs(thisdat.motionstrength);
        stim1(logical([(diff(thisdat.trialnr) ~= 1); 1])) = NaN;
        [rho, pval] = corr(stim1, circshift(stim1, 1), 'rows', 'complete');
        stimRep.rho(sj, session) = rho;
        stimRep.pval(sj, session) = pval;
    end
end
fprintf('rho = %.3f, range %.3f-%.3f, significant in %d out of %d sessions', ...
    nanmean(stimRep.rho(:)), min(stimRep.rho(:)), max(stimRep.rho(:)), ...
    length(find(stimRep.pval(:) < 0.05)), numel(stimRep.pval));

% 4. mediation analysis

%% there you go! get in touch if you have any further questions.

% Anne Urai, anne.urai@gmail.com / @AnneEUrai
