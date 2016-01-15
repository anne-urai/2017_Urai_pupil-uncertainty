
sigmoid = @(x, mu) 1 ./ (1 + exp(-x + mu));
x = -6:0.01:6;

subplot(661); % tiny!
hold on;
plot(x, sigmoid(x, 0), 'k');
plot(x, sigmoid(x, 1.5), 'k-.');
plot(x, sigmoid(x, -1.5), 'k-.');


% axes
plot([-5 5], [0.5 0.5], 'k');
plot([0 0], [0 1], 'k');

xlim([-6 6]);
axis off;

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/boundShift.pdf'));
