function psychFuncShift_Bias_byEvidence(whichmodulator, nbins, correctness)
% instead of splitting by previous trial pupil, do this 3x for current
% trial evidence

% even better: make the same plot as Fruend, figure 1c and put in the
% supplement

global mypath;

% plot both the effect of pupil on overall repetition bias and show that
% this is symmetrical for both previous choices
subjects = 1:27;
clear grandavg;
nbins = 3;

for sj = unique(subjects),
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data(find(data.sessionnr > 1), :);
    
    % code for choice repetition from the last trial to this one
    repetition          = (diff(data.resp) == 0);
    data.repeat         = [NaN; repetition];

    [grandavg.ev(sj, :), grandavg.repeat(sj, :)] = ...
        divideintobins(abs(data.motionstrength), data.repeat, nbins);
end

% compute the deviation from 50%


subplot(441);
boundedline(1:nbins, mean(grandavg.repeat), std(grandavg.repeat) ./ sqrt(length(subjects)));

% also try with absolute influence of history?
