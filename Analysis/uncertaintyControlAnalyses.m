%% reproduce all the possible control analyses with model-based uncertainty I can think of
clear
global mypath;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

% fixed effects
%data = addUncertainty2data;
data = readtable(sprintf('%s/Data/CSV/2ifc_data_allsj.csv', mypath));
data.evidence       = data.coherence;
repetition          = (diff(data.resp) == 0);
data.repeat         = [repetition; NaN];

CorrelationMatrixUncertainty(data);
% scatter3dBinned(data); % explore this in 3d

% random effects
% addUncertainty2CSV; % only needs to be done once
uncertaintyRepetition_grandaverage;
% uncertaintyRepetition_individual;

% fruend model - run this once
mods{1} = 'pupil+uncertainty';
cd(sprintf('%s/Code/serial-dependencies', mypath));
system([sprintf('for sj in {1..27}; do filename=$(printf "data/2ifc_%s_sj%%02d.txt" $sj);', mods{1}), ...
    sprintf('echo $filename; python2.7 analysis.py -fr -n10 -p "%s/Data/serialmodel/" $filename; sleep 5; done', mypath)]);
a6_retrieveDataFromPython(mods{1}); % outputs absolute errors

FruendModelUncertainty; % 

%% fit a model where each of the 3 measures of uncertainty is included
HistoryModel3Regressors;
