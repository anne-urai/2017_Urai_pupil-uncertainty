function err = cumWB_LL(p, intensity, responses)
% see http://courses.washington.edu/matlab1/Lesson_5.html#1

% compute the vector of responses for each level of intensity
w = Weibull(p, intensity);

% negative loglikelihood, to be minimised
err = -sum(responses .*log(w) + (1-responses).*log(1-w));

end