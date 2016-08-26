
% unbiased logistic
x = linspace(-6, 6, 100);
logFunc = 1 ./ (1 + exp(-x));
close;
subplot(221);
plot(x, logFunc, 'b');

% two separate ones that are biased towards left or right
logFuncLeft = 1 ./ (1 + exp(-x + 2));
logFuncRight= 1 ./ (1 + exp(-x - 2));

hold on;
plot(x, logFuncLeft, 'g', x, logFuncRight, 'c');
grid on;

% compare an unbiased function with one that is the average of two biased ones
subplot(222);
logFuncMean = mean([logFuncLeft; logFuncRight]);
plot(x, logFunc, 'b', x, logFuncMean, 'r');
grid on;