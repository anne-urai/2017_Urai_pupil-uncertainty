function PsychometricFunction

global mypath;
nbins = 6;
subjects = 1:27;

grandavg.accuracy = nan(length(subjects), nbins);
grandavg.rt       = nan(length(subjects), nbins);

for sj = subjects,
    
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % evidence strength
    data.motionstrength = (abs(data.motionstrength));
    
    % bin
    [~, grandavg.accuracy(sj, :)] = divideintobins(abs(data.motionstrength), data.correct, nbins);
    [~, grandavg.rt(sj, :)] = divideintobins(data.motionstrength, data.rt, nbins);

end

% plot psychometric func and chronometric func in 1
[ax, hLine1, hLine2] = plotyy(1:nbins, 100*squeeze(nanmean(grandavg.accuracy)), ...
    1:nbins, squeeze(nanmean(grandavg.rt)));
hLine1.LineStyle = 'none';
hLine2.LineStyle = 'none';

cols(1, :) = [0.2 0.2 0.2];
cols(2, :) = [0.6 0.6 0.6];

% add errorbars
hold(ax(1), 'on');
errorbar(ax(1), 1:nbins, 100*squeeze(nanmean(grandavg.accuracy)), ...
    100*squeeze(nanstd(grandavg.accuracy)) ./ sqrt(length(subjects)), ...
    'color', cols(1,:), 'linewidth', 1);

set(ax(1), 'xlim', [0.5 nbins+0.5], 'xtick', [1 3.5 6], 'ylim', [50 100], ...
    'ytick', 50:25:100, 'box', 'off', 'ycolor', cols(1,:), 'xticklabel', {'weak', 'medium', 'strong'}); 
xlabel(ax(1), 'Evidence');
ylabel(ax(1), 'Accuracy (%)');
axis(ax(1), 'square');

hold(ax(2), 'on');
errorbar(ax(2), 1:nbins, squeeze(nanmean(grandavg.rt)), ...
    squeeze(nanstd(grandavg.rt)) ./ sqrt(length(subjects)), 'color', cols(2,:), 'linewidth', 1);

set(ax(2), 'xlim', [0.5 nbins+0.5], 'xtick', [1 3.5 6], 'ylim', [0.25 0.55], ...
    'ytick', 0.25:0.15:0.55, 'box', 'off', 'ycolor', cols(2,:)); 
xlabel(ax(2), 'Evidence');
ylabel(ax(2), 'Reaction time (s)');
axis(ax(2), 'square');

ax(2).YLabel.Rotation = 270;
ypos = ax(2).YLabel.Position;
ypos(1) = ypos(1)+0.8;
ax(2).YLabel.Position = ypos;

axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes),
  axes(a).FontSize = 7;
end


end