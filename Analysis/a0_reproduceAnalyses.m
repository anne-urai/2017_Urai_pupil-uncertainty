% 

% first, make sure this path matches the place where the data are stored
global mypath;
mypath = '~/Data/pupilUncertainty';
cd(mypath);

% add the code folder
addpath('~/Dropbox/Code/pupilUncertainty/Analysis');

% optional: get my Tools folder here https://github.com/anne-urai/Tools and
% add this to your path as well - or, manually get those functions that
% you're missing (mainly for plotting and some stats). 

% install fieldtrip - see http://www.fieldtriptoolbox.org 
addpath('~/Documents/fieldtrip/');
ft_defaults;

%% process pupil data
for sj = 1:27, a1_PupilAnalysis(sj); end

%% motion energy filtering
% warning: since this is memory and cpu heavy, the outputs of this script
% have been provided in motionEnergy.zip - consider using those to continue
% if you don't have access to a computing cluster.
for sj = 1:27, a2_MotionEnergy(sj); end

%% combine all the data that was processed into one grand average file
a3_writeData2GA;

%% also write one CSV file per sj for quick behavioural data analysis
% note: these files are also included in CSV.zip, so you could skip
% everything up to here if you're not interested in visualizing the pupil
% timecourses 
a4_writeData2CSV;

%% reproduce figure 1, 2, and 3
fig1ab_Model;
fig2_overview;
fig3_overview;

%% run the python model with modulatory term
% write away text files that have the format Python needs
a5_writeFiles4pythonToolbox;

% this is easiest to run from the terminal. With Python 2.7 installed, go
% to the folder mypath/Code/serial-dependencies 
cd(sprintf('%s/Code/serial-dependencies', mypath));

% then, call the terminal from Matlab (not sure if this would work on Windows)
mods = {'plain', 'pupil-rt', 'rt-pupil', 'feedbackpupil', 'fb-decpupil'};
if ~exist(sprintf('%s/Data/serialmodel', mypath), 'dir'), mkdir(sprintf('%s/Data/serialmodel', mypath)); end

for m = 1:length(mods),
    system([sprintf('for sj in {1..27}; do filename=$(printf "data/2ifc_%s_sj%%02d.txt" $sj);', mods{m}), ...
       sprintf('echo $filename; python2.7 analysis.py -fr -n1000 -p "%s/Data/serialmodel" $filename; sleep 5; done', mypath)]);
end
% important: if you want the quick and dirty version without accurate errorbars, change -n1000 to n-10.

%% get the output back into something Matlab can work with
cd(sprintf('%s/Code/Analysis', mypath));
for m = 1:length(mods),
    a6_retrieveDataFromPython(whichmodulator);
end

%% reproduce figure 4
fig4_overview; 

%% also reproduce the supplementary figures
figS1_MotionEnergy_Filters;
figS2_decisionStrategies_pupilSplit;
figS3_feedbackpupil;

% there you go! get in touch if you have any further questions.
% Anne Urai, anne.urai@gmail.com
