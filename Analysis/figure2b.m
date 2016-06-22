

% new figure with more uncertainty signatures
close; figure;
global mypath;

subplot(441); ReactionTimeUncertainty;

% other metrics of uncertainty:
subplot(443); UncertaintyAccuracy('rt');
subplot(444); PsychFuncs_byUncertainty('rt');

print(gcf, '-dpdf', sprintf('%s/Figures/figure2b.pdf', mypath));
