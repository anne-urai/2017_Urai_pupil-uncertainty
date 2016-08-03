function [slope, threshold, lapse] = fitWeibull(x,y)

pBest = fminsearchbnd(@(p) cumWB_LL(p, ...
    x, y), ...
    [1 3 0.1], [0 0 0], [3 6 1]);

slope       = pBest(1);
threshold   = pBest(2);
lapse       = pBest(3);

end

function err = cumWB_LL(p, intensity, responses)
% see http://courses.washington.edu/matlab1/Lesson_5.html#1

% compute the vector of responses for each level of intensity
w   = Weibull(p, intensity);

% negative loglikelihood, to be minimised
err = -sum(responses .*log(w) + (1-responses).*log(1-w));

end