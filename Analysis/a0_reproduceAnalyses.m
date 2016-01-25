% 

% first, make sure this path matches the place where the data are stored
global mypath;
mypath = '~/Desktop/newPupilUncertainty';
cd mypath;

% add the code folder
addpath('/Code');

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

%% also already write away text files that have the format Python needs
a5_writeFiles4pythonToolbox;

%% let's have a look at the model
% pdfs will be saved into mypath/Figures
fig1ab_Model;

