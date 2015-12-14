function fig4g_SjCorrelation(lagGroups, whichmodulator)
if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

clc;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));


% ============================================ %
% combine lags
% ============================================ %
flds = fieldnames(dat);

for f = 1:length(flds),
    try
        groupedDat.(flds{f}) = mean(dat.(flds{f})(:, lagGroups), 2);
    end
end

plot(abs(groupedDat.response), groupedDat.response_pupil, '.');

[rho, pval ] = corr(abs(groupedDat.response),  groupedDat.response_pupil);
title({sprintf('rho = %.3f', rho); sprintf('p = %.3f', pval)});

%axis tight;
axis square; xlabel('| Response weight |'); ylabel('Pupil interaction');
box off;