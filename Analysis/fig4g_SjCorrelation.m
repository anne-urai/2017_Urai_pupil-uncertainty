function fig4g_SjCorrelation(lagGroups, whichmodulator)

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

clc;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));
if ~exist('lagGroups', 'var'), lagGroups = 1; end

subplot(4,4,12);

% ============================================ %
% combine lags
% ============================================ %
flds = fieldnames(dat);

for f = 1:length(flds),
    try
        groupedDat.(flds{f}) = mean(dat.(flds{f})(:, lagGroups), 2);
    end
end


plot(groupedDat.response, groupedDat.response_pupil, '.');
%axis tight;
axis square; xlabel('Response weights'); ylabel('Pupil interaction');
box off;