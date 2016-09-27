% This code reproduces the analyses in the paper
% Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven 
% by decision uncertainty and alters serial choice bias. 
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai, 2016
% anne.urai@gmail.com

%% make sure this path matches the place where the data are stored
% and you unzipped everything
clear all; close all; clc;
global mypath;

% determine the path to the data
usr = getenv('USER');
switch usr
    case 'aurai' % uke cluster
        mypath = '~/Data/pupilUncertainty';
        addpath('~/code/pupilUncertainty/Analysis/');
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
for sj = 1:27, a1a_PupilAnalysis_NaNs(sj); end
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

%% plot figures 1-4
cd(mypath); if ~exist('Figures', 'dir'); mkdir Figures; end
figure1;
figure2;
figure3; 
% note: to get Bayes Factors, need R installed!
figure4; 

%% run the python model with modulatory term
% write away text files that have the format Python needs
a5_writeFiles4pythonToolbox;

% then, call the terminal from Matlab (not sure if this would work on Windows)
mods = {'plain', 'plainCoh', 'pupil+rt', 'fbpupil', 'fb+decpupil', 'pupil', 'rt'};
if ~exist(sprintf('%s/Data/serialmodel', mypath), 'dir'), mkdir(sprintf('%s/Data/serialmodel', mypath)); end

% this is easiest to run from the terminal. With Python 2.7 installed, go
% to the folder mypath/Code/serial-dependencies
cd(sprintf('%s/Code/serial-dependencies', mypath));

for m = 1:length(mods),
    system([sprintf('for sj in {1..27}; do filename=$(printf "data/2ifc_%s_sj%%02d.txt" $sj);', mods{m}), ...
        sprintf('echo $filename; python2.7 analysis.py -fr -n1000 -p "%s/Data/serialmodel/" $filename; sleep 5; done', mypath)]);
end

% important: if you want the quick and dirty version without accurate
% individual errorbars, change -n1000 to -n10. Otherwise, the code will run 
% bootstraps which can take quite a while.

% get the output back into something Matlab can work with
cd(sprintf('%s/Code/Analysis', mypath));
for m = 1:length(mods),
    a6_retrieveDataFromPython(mods{m}); 
end

close all;
figure5;

%% also create the supplementary figures
figureS1;
figureS2;
figureS3;
figureS4;
figureS5;
figureS6;
figureS7;
figureS8;
figureS9;
figureS10; 
figureS11;

%% analyze some final stuff that's not in the figures

% reported in text
stimulusTransitions;

% do individual choice tendency correlate with psychometric thresholds?
data          = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
data.evidence = abs(data.motionstrength);
pBest         = rowfun(@fitWeibull, data, 'inputvariables', {'evidence', 'correct'}, ...
    'groupingvariables', {'subjnr'}, 'outputvariablenames', {'slope', 'threshold', 'lapse'});
load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, 'plain'));
scatter(pBest.threshold', dat.response(:, 1));
[rho, pval] = corr(pBest.threshold, dat.response(:, 1));

% RT for the easiest trials
data    = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
data    = data(find(data.difficulty == 5), :);
data.rt = data.rt * 1000;
RTs     = splitapply(@median, data.rt, findgroups(data.subjnr)); % in ms
fprintf('RTs for the easiest trials: mean %.3f, min %.3f, max %.3f', mean(RTs), min(RTs), max(RTs));

% does figure 5 depend on running pupil and RT in the same model?
subplot(5,5,1); rho1 = SjCorrelation('pupil', 'response', 'pupil');
ylabel('Pupil x choice weight');
subplot(5,5,2); rho2 = SjCorrelation('pupil', 'response', 'rt');
ylabel('RT x choice weight'); set(gca, 'yaxislocation', 'right');
suplabel('Choice weight', 'x');
[ridiff,cilohi,p] = ridiffci(rho1, rho2, 27, 27, 0.05);
suplabel(sprintf('delta r = %.3f, p = %.3f', ridiff, p), 't');

% median of ISI
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj_withlatencies.csv', mypath));
ISIs = splitapply(@nanmedian, data.latency_feedbackisi, findgroups(data.subjnr)) - 2.25;

% does the pupil interpolation affect the estimated pupil x choice weight?
load(sprintf('%s/Data/GrandAverage/historyweights_pupil+rt.mat', mypath));
load(sprintf('%s/Data/GrandAverage/pupilquality.mat', mypath));
[rho, pval] = corr(dat.response_pupil(:, 1), grandavg.rejectTrials');

%% there you go! get in touch if you have any further questions.
% Anne Urai, anne.urai@gmail.com / @AnneEUrai
