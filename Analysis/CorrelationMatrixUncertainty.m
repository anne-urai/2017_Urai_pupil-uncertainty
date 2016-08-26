function CorrelationMatrixUncertainty(data)
global mypath;

% =================================================== %
% correlate all to all!
% =================================================== %
close all;
nbins = 3;
whichVars = {'evidence', 'decision_pupil', 'rt', 'repeat'};
%whichVars = {'evidence', 'decision_pupil', 'repeat'};

colors = [cbrewer('qual', 'Set1', 2); [0 0 0]];
colors = colors([1 2 4], :);

for correctness = [0 1 2],
    
    if correctness < 2,
        thisdat   = data((data.correct == correctness), :);
    else
        thisdat = data;
    end
    nsubpl    = length(whichVars);
    
    for i = 1:nsubpl,
        for j = 1:nsubpl,
            
            % start a new subplot
            subplot(nsubpl, nsubpl, i + (j-1)*nsubpl);
            hold on; % for later correct and error
            
            if i == j,
                % for autocorrelation, plot histfit
                h = histogram(thisdat.(whichVars{i}), ...
                    'edgecolor', 'none', 'facecolor', colors(correctness+1, :));
                axis tight; box off; axis square;
                
            elseif i < j,
                
                % bin the data rather than a full correlation!
                [binnedx, binnedy, stdx, stdy] = ...
                    divideintobins(thisdat.(whichVars{i}), thisdat.(whichVars{j}), nbins);
                h = ploterr(binnedx, binnedy, [], [], ...
                    '-o','hhxy',0.1);
                set(h(1), 'MarkerSize', 3, ...
                    'MarkerEdgeColor', colors(correctness+1, :), ...
                    'MarkerFaceColor', 'w', 'color', colors(correctness+1, :));
                
                % make the axes a little bit larger
                axis tight; axis square; axisNotSoTight;
                
                % get correlations
                if strcmp(whichVars{j}, 'repeat') && j ~= i,
                    % logistic beta of x variable onto repetition
                    mdlSpec = sprintf('repeat ~ 1 + %s', whichVars{i});
                    mdl = fitglm(thisdat, mdlSpec, 'distr', 'binomial', 'link', 'logit');
                    rho(i + (j-1)*nsubpl, correctness+1) = mdl.Coefficients.Estimate(2);
                else
                    % linear beta
                    mdlSpec = sprintf('%s ~ 1 + %s', whichVars{j}, whichVars{i});
                    mdl = fitglm(thisdat, mdlSpec);
                    rho(i + (j-1)*nsubpl, correctness+1) = mdl.Coefficients.Estimate(2);
                end
                if correctness == 1,
                    title(['{\color{red}\beta ' sprintf('%+.2f', rho(i + (j-1)*nsubpl, 1)) '} {\color{blue}\beta ' ...
                        sprintf('%+.2f', rho(i + (j-1)*nsubpl, 2)) '}']);
                end
            else
                % leave white, dont plot the correlations twice
                axis off;
            end
            
            % layout
            if j == nsubpl,               xlabel(whichVars{i}, 'interpreter', 'none'); end
            if i == 1,                    ylabel(whichVars{j}, 'interpreter', 'none'); end
            if j < nsubpl && j ~= i,      set(gca, 'xtick', []); end
            if i > 1 && j ~= i,           set(gca, 'ytick', []); end
            
            if strcmp(whichVars{j}, 'repeat') && j ~= i, ...
                    set(gca, 'ylim', [0.47 0.56], 'ytick', [0.48:0.02:0.56], 'ygrid', 'on');
            end
            set(gca, 'tickdir', 'out', 'box', 'off');
        end
    end
end
set(findall(gcf,'type','text'),'FontSize',10);
print(gcf, '-dpdf', sprintf('%s/Figures/uncertaintyCorrMatrix_bins%d.pdf', mypath, nbins));

end