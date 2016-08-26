function out = ranktransform(in)
% takes a vector and returns a rank-transformed version

% in  = normrnd(0,1,10,1);
% idx = 1:length(in);

% default is from small to large
[sorted, sortidx] = sort(in);

% sortedidx now indicates in which order in needs to be to be sorted
% however, we want to know which place in the rank order each element of in
% needs to have!

out = 1:length(in);
out(sortidx) = out;

end