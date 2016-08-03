function unc = simulateUncertainty(x, correct, sigma, bound, n)
% compute uncertainty for correct and error trials based on stimulus strength
% x should be absolute, so that only one 'side' is simulated

if ~exist('n', 'var'), n = 10000; end

% simulate decision variables for this level of stimulus strength
dvs = x - bound + sigma * randn(n, 1);

if correct == 1,
    % find the mean distance to bound for the correct samples
    unc = mean(1 - tanh(abs(dvs(dvs > bound) - bound)));
elseif correct == 0,
    % find the mean distance to bound for error samples
    unc = mean(1 - tanh(abs(dvs(dvs < bound) - bound)));
end

end