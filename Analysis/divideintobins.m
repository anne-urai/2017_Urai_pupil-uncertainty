function [binnedx, binnedy, stdx, stdy, rho, pval] = divideintobins(x, y, nbins, corrtype, summaryFunc)
% This code reproduces the analyses in the paper
% Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven 
% by decision uncertainty and alters serial choice bias. 
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai, 2016
% anne.urai@gmail.com

% take two vectors and divide x them into nbins, then compute the mean
% response of y for each x

if ~exist('corrtype', 'var'); corrtype = 'Spearman'; end
if ~exist('nbins', 'var'); nbins = 3; end

%assert(~any(isnan(x)), 'x contains nans');
%assert(~any(isnan(y)), 'y contains nans');
if ~exist('summaryFunc', 'var'); summaryFunc = @nanmean; end % to allow for median

switch func2str(summaryFunc)
    case 'nanmean'
        distFun = @nanstd;
    case {'nanmedian', 'median'}
        distFun = @iqr;
end

if nbins == length(unique(x)),
    % split into the categories of x
    thisx = unique(x);
    mybins = nan(1, length(thisx) - 1);
    for i = 1:length(thisx) -1,
        mybins(i) = mean(thisx(i:i+1));
    end
    binIdx = discretize(x, [-inf mybins inf]);
elseif nbins == 2,
    binIdx = discretize(x, [-inf median(x) inf]);
else
    binIdx = discretize(x, [-inf quantile(x, nbins-1) inf]);
end

% get the summary measure for each bin
binnedx = splitapply(summaryFunc, x, binIdx);
binnedy = splitapply(summaryFunc, y, binIdx);

% also get the distribution within each bin
stdx = splitapply(distFun, x, binIdx);
stdy = splitapply(distFun, y, binIdx);

% do some statistics
if ~isempty(corrtype),
    % on binned or non-binned data? binning will change rho, but not really 
    % the p-value of the correlation
    [rho, pval] = corr(x(:), y(:), 'type', corrtype);
else
    rho = []; pval = [];
end

end