function y = Weibull(p, x)
% Parameters:  p(1) slope
%             p(2) threshold yeilding ~80% correct
%             p(3) lapse rate
%             x   intensity values.

g = 0.5;  %chance performance

% include a lapse rate, see Wichmann and Hill parameterisation
y = g + (1 - g - p(3)) * (1-exp(- (x/p(2)).^p(1)));

end