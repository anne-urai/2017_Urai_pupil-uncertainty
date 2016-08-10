% Urai AE, Braun A, Donner TH. Decision uncertainty drives pupil-linked
% arousal systems and alters subsequent decision-making.

% first, make sure this path matches the place where the data are stored
% and you unzipped everything
global mypath;
% determine the path to the data
usr = getenv('USER');
switch usr
    case 'aurai' % uke cluster
        mypath = '~/Data/pupilUncertainty';
    case 'anne' % macbook pro
        mypath = '/Users/anne/Data/pupilUncertainty_FigShare';
        addpath(genpath('~/Dropbox/code/Tools/'));
end

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

% check how much of the data is interpolated
for sj = 16:27, a1a_PupilAnalysis_NaNs(sj); end
pupilInterpolationCount; % compute a percentage per trial and make an overview

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
mods = {'plain', 'pupil+rt', 'fbpupil', 'fb+decpupil', 'pupil', 'rt'};
if ~exist(sprintf('%s/Data/serialmodel', mypath), 'dir'), mkdir(sprintf('%s/Data/serialmodel', mypath)); end

% this is easiest to run from the terminal. With Python 2.7 installed, go
% to the folder mypath/Code/serial-dependencies
cd(sprintf('%s/Code/serial-dependencies', mypath));

for m = 1:length(mods),
    system([sprintf('for sj in {1..27}; do filename=$(printf "data/2ifc_%s_sj%%02d.txt" $sj);', mods{m}), ...
        sprintf('echo $filename; python2.7 analysis.py -fr -n10 -p "%s/Data/serialmodel/" $filename; sleep 5; done', mypath)]);
end

% important: if you want the quick and dirty version without accurate
% individual errorbars, change -n1000 to -n10. Otherwise, the code will run 
% bootstraps which can take quite a while.

%% get the output back into something Matlab can work with
cd(sprintf('%s/Code/Analysis', mypath));
for m = 1:length(mods),
    a6_retrieveDataFromPython(mods{m}); 
end

%% reproduce figures 4,5,6 based on the model fits
close all;
figure4;

close all;
figure5;

% confirm that pupil only changes sequential biases
postErrorSlowing;

close all; 
figure6;

%% also create the supplementary figures
figureS2_RT;
figureS4_History_CorrectError;
figureS3_feedbackpupil;
figureS5_pupilResponseLagged;
figureS5_scatterIndividual; 
figureS6_performanceOverSessions;
mediationAnalysis;
figureS1_MotionEnergy_Filters;

%% analyze some final stuff that's not in the figures
% reported in text, or only shown in the response to reviewers

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

% 3. compute autocorrelation in evidence strength
stimRep.rho = nan(27, 7);
stimRep.pval = nan(27, 7);
for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    for session = unique(data.sessionnr)',
        thisdat     = data(find(session==data.sessionnr), :);
        stim1       = abs(thisdat.motionstrength);
        stim1(logical([(diff(thisdat.trialnr) ~= 1); 1])) = NaN;
        [rho, pval] = corr(stim1, circshift(stim1, 1), 'rows', 'complete');
        stimRep.rho(sj, session) = rho;
        stimRep.pval(sj, session) = pval;
    end
end
fprintf('rho = %.3f, range %.3f-%.3f, significant in %d out of %d sessions', ...
    nanmean(stimRep.rho(:)), min(stimRep.rho(:)), max(stimRep.rho(:)), ...
    length(find(stimRep.pval(:) < 0.05)), numel(stimRep.pval));

% 5. checks on model-based uncertainty and its effect on switching
uncertaintyControlAnalyses;

% 6. does individual choice tendency correlate with psychometric thresholds?
data          = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
data.evidence = abs(data.motionstrength);
pBest         = rowfun(@fitWeibull, data, 'inputvariables', {'evidence', 'correct'}, ...
    'groupingvariables', {'subjnr', 'sessionnr'}, 'outputvariablenames', {'slope', 'threshold', 'lapse'});

load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));
scatter(pBest(:, 2), dat.response(:, 1));
[rho, pval] = corr(pBest(:, 2), dat.response(:, 1));

% 7. RT for the easiest trials
data    = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
data    = data(find(data.difficulty == 5), :);
data.rt = data.rt * 1000;
RTs     = splitapply(@median, data.rt, findgroups(data.subjnr)); % in ms
fprintf('mean %.3f, min %.3f, max %.3f', mean(RTs), min(RTs), max(RTs));

% 8. baseline pupil does not predict repetition behaviour
clf; subplot(441); b = Uncertainty_byErrorCorrect('baseline_pupil');
subplot(4,7,3); plotBetasSwarm(b, colors([1 2], :));
set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'});

subplot(443); psychFuncShift_Bias('baseline_pupil', 3, []);
print(gcf, '-dpdf', sprintf('%s/Figures/baselinePupil.pdf', mypath));

% 9. does figure 6 depend on running pupil and RT in the same model?
subplot(5,5,1); rho1 = SjCorrelation('pupil', 'response', 'pupil');
ylabel('Pupil * choice weight');
subplot(5,5,2); rho2 = SjCorrelation('pupil', 'response', 'rt');
ylabel('RT * choice weight'); set(gca, 'yaxislocation', 'right');
suplabel('Choice weight', 'x');

% unpaired test between rho's (because the fitted choice weights are not
% identical anymore)
[ridiff,cilohi,p] = ridiffci(rho1, rho2, 27, 27, 0.05);
suplabel(sprintf('delta r = %.3f, p = %.3f', ridiff, p), 't');
print(gcf, '-dpdf', sprintf('%s/Figures/separateRegressionModels.pdf', mypath));

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

subplot(4,4,[1 2]); hold on;
colors = cbrewer('seq', 'Greens', 10); colors = colors(8:end, :);
boundedline(dat.newvariance.stimuli, ...
    squeeze(nanmean(dat.newvariance.explained)), ...
    squeeze(nanstd(dat.newvariance.explained)) ./ sqrt(sum(~isnan(dat.newvariance.explained))), ...
    'cmap', colors(1, :));
xlabel('\Delta motion coherence'); ylabel('History contribution (%)');

print(gcf, '-dpdf', sprintf('%s/Figures/historyContribution.pdf', mypath));

%% there you go! get in touch if you have any further questions.
% Anne Urai, anne.urai@gmail.com / @AnneEUrai
