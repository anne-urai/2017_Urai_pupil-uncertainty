function boundUpdate_Kepecs
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
% Reproduces the heatmap figure from Kepecs 2007 COSYNE poster
%
% Anne Urai, 2017
% anne.urai@gmail.com

clearvars -except mypath; close all;
global mypath;
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

% warnings
warning('error', 'stats:glmfit:PerfectSeparation');
warning('error', 'stats:glmfit:IterationLimit');
warning('error', 'stats:glmfit:IllConditioned');

% get all data
prevbins    = 8; % has to be even
colormap(cbrewer('div', 'RdBu', 64));
set(groot, 'defaultaxescolororder', viridis(3), 'DefaultAxesFontSize', 7);

cors = [1 0];
for c = 1:2,
    clf;
    
    subjects = 1:27;
    for sj = subjects,
        
        % get data
        data     = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        % make ME evenly distributed over trials
        for s = unique(data.sessionnr)',
            data.motionstrength(s == data.sessionnr) = ...
                zscore(data.motionstrength(s == data.sessionnr));
        end
        
        % define by motion strength
        data.prevcorrect        = circshift(data.correct, 1);
        data.prevmotionstrength = circshift(data.motionstrength, 1);
        data.prevmotionbin      = discretize(data.prevmotionstrength, ...
            [-inf quantile(data.prevmotionstrength, prevbins-1) inf]);

        % define by difficulty
        % data.prevmotionbin      = circshift(discretize(data.difficulty .* data.stim, [-6 -3.5 -1.5 0 1.5 3.5 6]), 1);
        % data.difficulty         = discretize(abs(data.difficulty), [0 2 6]);
        qt = quantile(abs(data.motionstrength), 2); % make s
        data.currDifficulty       = discretize(abs(data.motionstrength), [0 qt(1) Inf]);

        % ==================================================================
        % split and fit bias
        % ==================================================================
        
        [gr, currIdx, prevDiff, prevCorr] = findgroups(data.currDifficulty, ...
            data.prevmotionbin, data.prevcorrect);
        
        bias = splitapply(@fun, data.motionstrength, (data.resp > 0), gr);
        bias = cat(2, bias{:});
        
        intercept           = bias(1, :)'; % get bias term
        slope               = bias(2, :)'; % get bias term
        
        dat{sj} = array2table([currIdx prevDiff prevCorr sj*ones(size(currIdx)) intercept slope], ...
            'variablenames', {'currEv', 'prevEv', 'prevCorr', 'subj_idx', 'bias', 'slope'});
        
        % fill with NaNs if some are missing
        try
            assert(numel(unique(gr)) == prevbins*2*2);
        catch
            disp('adding nans');
            % which combinations should there be?
            sets        = {unique(currIdx), unique(prevDiff), unique(prevCorr)};
            [x, y, z]   = ndgrid(sets{:});
            cartProd    = [x(:) y(:) z(:)];
            missing     = find(~ismember(cartProd, dat{sj}{:, 1:3}, 'rows'));
            addStuff      = array2table([cartProd(missing, 1) cartProd(missing, 2) cartProd(missing, 3) ...
                sj*ones(size(missing)) nan(size(missing)) nan(size(missing))], ...
                'variablenames', {'currEv', 'prevEv', 'prevCorr', 'subj_idx', 'bias', 'slope'});
            dat{sj}     = [dat{sj}; addStuff];
        end
        
        subplot(5,6,sj); hold on;
        newx = linspace(-1, 1, numel(unique(prevDiff)));
        plot(newx, zeros(size(newx)), 'color', [0.5 0.5 0.5]);
        colors = cbrewer('qual', 'Set2', 2);
        
        for currEv = 1:2,
            thisdat     = dat{sj}.bias( dat{sj}.currEv == currEv & dat{sj}.prevCorr == cors(c));
            plot(newx, thisdat, 'color', colors(currEv, :));
            plot(newx, thisdat, '.', 'color', colors(currEv, :));
            
        end
        title(sprintf('P%02d', sj)); box off;
        ylims = get(gca, 'ylim');
        ylims = max(abs(ylims));
        ylim([-ylims ylims]);
    end
    
    dat2     = cat(1, dat{:});
    dat2.bias(abs(dat2.bias) > prctile(abs(dat2.bias), 97)) = NaN;
    % dat2    = unstack(dat, {'bias', 'slope'}, {'subj_idx'});
    
    subplot(5,6,28); hold on;
    plot(newx, zeros(size(newx)), 'color', [0.5 0.5 0.5]);
    
    prevbins = numel(unique(prevDiff));
    for currEv = 1:2,
        thisdat     = reshape(dat2.bias(dat2.currEv == currEv & dat2.prevCorr == cors(c)), prevbins, 27)';
        p(currEv)   = boundedline(newx, squeeze(nanmean(thisdat)), ...
            squeeze(nanstd(thisdat)) ./ sqrt(numel(subjects)), 'alpha', 'cmap', colors(currEv, :));
        hold on;
        plot(newx, squeeze(nanmean(thisdat)), '.', 'color', colors(currEv, :), 'markersize', 5);
    end
    title('Grand Average');
    box off;
    
    % if c == 2,
    l = legend(p, {'Current trial hard', 'Current trial easy'});
    legend boxoff;
    l.Position(1) = l.Position(1) + 0.2;
    
    switch cors(c)
        case 0
            ylim([-0.5 0.5]);
            suplabel('Previous error', 't');
            name = 'error';
        case 1
            ylim([-0.5 0.5]);
            suplabel('Previous correct', 't');
            name = 'correct';
    end
    
    suplabel('Previous stimulus strength', 'x');
    suplabel('Choice bias', 'y');
    
    print(gcf, '-dpdf', sprintf('%s/Figures/boundUpdate_%s.pdf', mypath, name));
end
end


function b = fun(x1, x2)

% skip if any errors were thrown
try
    b = glmfit(x1, x2, 'binomial', 'link', 'logit');
    
    % include lapse rates
    %[b1, b2] = fitLogistic(x1, x2);
   % b = [b1; b2];
catch
    b = [NaN; NaN];
end

b = {b};
end


function [bias, slope, lapseLow, lapseHigh] = fitLogistic(x,y)

% make gamma and lambda symmetrical
pBest = fminsearchbnd(@(p) logistic_LL(p, ...
    x, y), [0 1 0.01], [-4 0 0], [4 10 0.1]);

bias        = pBest(1);
slope       = pBest(2);
lapseLow    = pBest(3);
lapseHigh   = pBest(3);

end

function err = logistic_LL(p, intensity, responses)
% see http://courses.washington.edu/matlab1/Lesson_5.html#1

% compute the vector of responses for each level of intensity
w   = logistic(p, intensity);

% negative loglikelihood, to be minimised
err = -sum(responses .*log(w) + (1-responses).*log(1-w));

end


