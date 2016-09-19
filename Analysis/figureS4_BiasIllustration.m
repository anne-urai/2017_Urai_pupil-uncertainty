
% supplementary figure
global mypath;
close; figure;

% first, illustrate shift
sigmoid = @(x, mu) 1 ./ (1 + exp(-x + mu));
x = -5:0.01:5;
subplot(441); % tiny!
hold on;
plot(x, sigmoid(x, 0), 'k');
plot([-5 5], [0.5 0.5], 'k');
plot([0 0], [0 1], 'k');
axis tight; axis square;
ylabel('P(choice = 1)');
xlabel('Sensory evidence (a.u.)');
ylim([-0.05 1]); 

subplot(442); ExamplePsychFunc(20); title(sprintf('Participant %d', 20));
subplot(443); ExamplePsychFunc(13); title(sprintf('Participant %d', 13));
subplot(444); ExamplePsychFuncShift(10); title('Participant 10');

print(gcf, '-dpdf', sprintf('%s/Figures/figureS4.pdf', mypath));
