
%

clear; close all; clc;

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

%% 
figure;
% total average and average split by accuracy x difficulty
fig3a_PupilTimecourse(1);

% uncertainty timecourses
subplot(4,4,9);
fig3b_PupilUncertaintyTimecourse(0);

% bias timecourses
subplot(4,4,13);
fig3f_pupilDprimeCriterionTimecourse()

%% save fig
print(gcf, '-dpdf', sprintf('%s/Figures/figure_cosyne.pdf', mypath));





