
% reproduces 
global mypath;

close;
figure;

% take P10 as an example
subplot(441); ExamplePsychFuncShift(10);

subplot(445); FruendKernels('plain');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

subplot(446); decisionStrategies('plain');
set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);

print(gcf, '-dpdf', sprintf('%s/Figures/figure5.pdf', mypath));
