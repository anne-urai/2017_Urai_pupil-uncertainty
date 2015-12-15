function fig4S_FruendPlainVsMod(lagGroups, whichmodulator)
if ~exist('lagGroups', 'var'), lagGroups = 1; end
if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end

% ============================================ %
% bargraphs for previous response and response * pupil regressors
% ============================================ %

clc;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'plain'));
plainDat = dat;
load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', whichmodulator));
modDat = dat;

% ============================================ %
% combine lags
% ============================================ %
flds = fieldnames(plainDat);

for f = 1:length(flds),
    try
        plainDatGrouped.(flds{f}) = mean(plainDat.(flds{f})(:, lagGroups), 2);
    end
end

% ============================================ %
% combine lags
% ============================================ %
flds = fieldnames(modDat);

for f = 1:length(flds),
    try
        modDatGrouped.(flds{f}) = mean(modDat.(flds{f})(:, lagGroups), 2);
    end
end

hold on;
plot(plainDatGrouped.response, modDatGrouped.response, '.');
axis tight; 
ylims = get(gca, 'ylim'); xlims = get(gca, 'xlim');
newlims = max([ylims xlims])*1.3;

plot([-newlims newlims], [-newlims newlims], 'k');
plot(plainDatGrouped.response, modDatGrouped.response, 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize',3);
ylim([-newlims newlims]); xlim([-newlims newlims]);
xlabel('Without modulation'); ylabel('With modulation'); 
axis square; 
