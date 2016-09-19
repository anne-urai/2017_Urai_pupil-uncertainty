
% stimulus proportions from Kepecs
a = linspace(99, 1, 1000000);
b = fliplr(a);
% c = linspace(-7, 7, 100000);
% d = log(a ./ (1 - a));

stimsRatio       = a ./ b;
stimsLogRatio    = log(stimsRatio);
stimPerc         = a * 100;

t = [a' b' (a>b)' stimsRatio' stimsLogRatio'];
subplot(331);
plot(t); legend('a', 'b', 'choice', 'ratio', 'logratio');
subplot(332);
plot([(a>b)' log(stimsRatio)']);
legend('choice', 'logratio');

% get uncertainty values
model.evs           = stimsLogRatio;
model.dvs           = model.evs + normrnd(0, 0.75, size(model.evs));
model.choice        = sign(model.dvs);
model.correct       = (model.choice == sign(model.evs));
model.confidence    = tanh(abs(model.dvs));
model.uncertainty   = 1 - model.confidence;

% plot
[unc, acc] = divideintobins(model.uncertainty, model.correct, 100);
hold on;
subplot(333);
plot(unc, 100*acc); xlabel('Uncertainty'); ylabel('Accuracy');
axis square; box off; xlim([0 1]); ylim([50 100]);

