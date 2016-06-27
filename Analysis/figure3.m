

% new figure with more uncertainty signatures
close; figure;
global mypath;

% use nice shades of red and green
colors = cbrewer('qual', 'Set1', 9);

% error vs correct
subplot(441); b = Uncertainty_byErrorCorrect('rt');
subplot(473); plotBetasSwarm(b(:, :, 2), colors([1 2], :));
set(gca, 'xtick', [1 2], 'xticklabel', {'Error', 'Correct'});

% other metrics of uncertainty
subplot(4,4,5); b = UncertaintyAccuracy('rt');

% show betas
subplot(4,7,10); plotBetasSwarm(b(:, 2), [0 0 0]);
set(gca, 'xtick', 1, 'xticklabel', []);
ylim([-0.6 0]);

% psychometric functions
subplot(4,4,9); b = PsychFuncs_byUncertainty('rt');
subplot(4,7,17); plotBetasSwarm(b(:, :, 2), [0.7 0.7 0.7; 0.2 0.2 0.2]);
set(gca, 'xtick', [1 2], 'xticklabel', {'fast', 'slow'});
xlabel('Reaction time'); ylabel('Threshold');
ylim([0 7]);

% ensure same axes proportions
ax = findobj(gcf, 'type', 'axes');
for a = 1:length(ax),
    pos = get(ax(a), 'position'); pos(4) = 0.15; set(ax(a), 'position', pos);
end

print(gcf, '-dpdf', sprintf('%s/Figures/figure3.pdf', mypath));
