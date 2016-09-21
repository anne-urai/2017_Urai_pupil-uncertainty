
close; figure;
global mypath;

%% add psychometric functions, chronometric functions
for s = 2:6,
    subplot(5,5,s-1);
    psychometricFunction(s);
end

if ~exist(sprintf('%s/Data/GrandAverage/historyweights_plain_session6.mat', mypath), 'file'),
    
    %% run history model for every session
    cd(sprintf('%s/Code/serial-dependencies/data', mypath));
    subjects = 1:27;
    
    % for supplement, do this per session
    for sj = subjects,
        disp(sj);
        
        for session = 2:6,
            data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
            
            if sj == 15,
                % this person did half a session 3 and half a session 6 + full session 7
                % change this so that each chunk has 10 blocks in it
                data.sessionnr(466:937) = 3;
                data.sessionnr(938:1425) = 4;
                data.sessionnr(1426:1921) = 5;
                data.sessionnr(1922:end) = 6;
            end
            
            data = data(find(data.sessionnr == session), :);
            
            % generate block nrs, NOT identical to session nrs! History effects
            % should not continue beyond a block
            blockchange = find(diff(data.trialnr) < 0);
            blocknrs = zeros(height(data), 1);
            for b = 1:length(blockchange)-1,
                blocknrs(blockchange(b)+1:blockchange(b+1)) = blocknrs(blockchange(b))+1;
            end
            blocknrs(blockchange(end)+1:end) = blocknrs(blockchange(end))+1;
            
            % no modulation, just history
            newdat = [ blocknrs data.sessionnr abs(data.motionstrength) (data.motionstrength > 0) (data.resp > 0)];
            
            dlmwrite(sprintf('2ifc_plain_session%d_sj%02d.txt', session, sj), ...
                newdat,'delimiter','\t','precision',4);
        end
    end
    
    %% run the actual model
    
    cd(sprintf('%s/Code/serial-dependencies', mypath));
    for session = 2:6,
        system([sprintf('for sj in {1..27}; do filename=$(printf "data/2ifc_plain_session%d_sj%%02d.txt" $sj);', session), ...
            sprintf('echo $filename; python2.7 analysis.py -fr -n10 -p "%s/Data/serialmodel" $filename; sleep 5; done', mypath)]);
    end
    
    % retrieve those in matlab
    cd(sprintf('%s/Code/Analysis', mypath));
    for session = 2:6,
        a6_retrieveDataFromPython(sprintf('plain_session%d', session)); % outputs absolute errors
    end
end

% plot the history kernels, only mean +- sem
for s = 2:6,
    subplot(5,5,s-1+10);
    load(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, sprintf('plain_session%d', s)));
    boundedline(1:7, mean(dat.response), std(dat.response) ./ sqrt(27), 'k');
    
    set(gca, 'ycolor', 'k', 'xtick', 1:7, 'ytick', [-0.2 0 0.2], ...
        'ylim', [-.2 .2], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
    ylabel('');
    if s == 2, ylabel('Choice weight'); else set(gca, 'yticklabel', []); end
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    xlabel('Lags'); axis square;
end

% plot history kernels, only individual data
for s = 2:6,
    subplot(5,5,s-1+15);
    fruendKernels(sprintf('plain_session%d', s), 'response');
    set(gca, 'ycolor', 'k', 'xtick', 1:7, 'ytick', [-0.7 0 0.7], ...
        'ylim', [-.75 .7], 'xlim', [0.5 7.5], 'box', 'off', 'xminortick', 'off');
    ylabel('');
    if s == 2, ylabel('Choice weight'); else set(gca, 'yticklabel', []); end
    set(gca, 'xcolor', 'k', 'ycolor', 'k', 'linewidth', 0.5);
    axis square;
end


print(gcf, '-dpdf', sprintf('%s/Figures/FigureS2.pdf', mypath));

% ====================================================================================== %
%% do repetition bias stats for each session, big ANOVA
% ====================================================================================== %

subjects = 1:27;
clear grandavg;
lag = 1; correctness = []; nbins = 3;

for sj = unique(subjects),
    for session = 2:6,
        data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
        
        data = data((data.sessionnr == session), :); % get rid of residual learning effects
        data.decision_pupil = projectout(data.decision_pupil, zscore(log(data.rt + 0.1)));
        % outcome vector need to be 0 1 for logistic regression
        data.resp(data.resp == -1) = 0;
        
        % get an overall logistic fit for normalization
        grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, ...
            data.resp, 'binomial','link','logit');
        
        % previous response
        resps = [0 1];
        for r = 1:2,
            
            % split into quantiles for each response
            if isempty(correctness),
                uncQs = quantile(data.decision_pupil(data.resp == resps(r)), nbins - 1);
            else
                uncQs = quantile(data.decision_pupil(data.resp == resps(r) & data.correct == correctness), nbins - 1);
            end
            
            % uncertainty bins
            for u = 1:nbins,
                
                % find the trials in this bin
                if isempty(correctness), % take all trials
                    switch u
                        case 1
                            trls = find(data.resp == resps(r) & data.decision_pupil <= uncQs(u));
                        case nbins % last one
                            trls = find(data.resp == resps(r) & data.decision_pupil > uncQs(u-1));
                        otherwise
                            trls = find(data.resp == resps(r) & ...
                                data.decision_pupil > uncQs(u-1) & data.decision_pupil <= uncQs(u));
                    end
                else % either correct or error trials
                    switch u
                        case 1
                            trls = find(data.resp == resps(r) & data.correct == correctness & data.(whichMod) <= uncQs(u));
                        case nbins % last one
                            trls = find(data.resp == resps(r) & data.correct == correctness & data.(whichMod) > uncQs(u-1));
                        otherwise
                            trls = find(data.resp == resps(r) & data.correct == correctness & ...
                                data.(whichMod) > uncQs(u-1) & data.(whichMod) <= uncQs(u));
                    end
                end
                
                % with this selection, take the trials after that
                laggedtrls = trls+lag;
                
                % exclude trial at the end of the dataset
                if any(laggedtrls > size(data, 1)),
                    trls(laggedtrls > size(data, 1)) = [];
                    laggedtrls(laggedtrls > size(data, 1)) = [];
                end
                
                % remove trials that dont match in block nr
                removeTrls = find(data.blocknr(laggedtrls) ~= data.blocknr(trls));
                laggedtrls(removeTrls) = [];
                
                % fit logistic regression
                thisdat = data(laggedtrls, :);
                
                try
                    b = glmfit(thisdat.motionstrength, thisdat.resp, ...
                        'binomial','link','logit');
                catch
                    warning('putting nan in betas!');
                    b = nan(size(b));
                end
                
                % save betas
                grandavg.logistic(sj, session-1, r, u, :) = b;
                
            end % uncertainty bin
        end % resp
    end % sj
end

% ========================================================= %
% normalize by their overall bias for A or B
% keep only the bias, discard slope info
% ========================================================= %

grandavg.logistic = grandavg.logistic(:, :, :, :, 1);

% ========================================================= %
% compute repetition (rather than response) bias
% ========================================================= %

resp1 = -1 * squeeze(grandavg.logistic(:, :, 1, :));
resp2 = squeeze(grandavg.logistic(:, :, 2, :));

% since this is centred at 0.5, treat it that way
grandavg.logisticRep = (resp1 + resp2) ./ 2;

% ========================================================= %
% transform from log(odds) to probability
% ========================================================= %

grandavg.logisticRep = exp(grandavg.logisticRep)./(1+exp(grandavg.logisticRep));

% ========================================================= %
% stats
% ========================================================= %

y1      = squeeze(grandavg.logisticRep);
s       = repmat(transpose(1:27), [1 5 3]);
ft      = repmat(1:5, [27 1 3]); % session
f1      = ft(:);
ft      = permute(repmat(1:3, [27 1 5]), [1 3 2]); % pupil
f2      = ft(:);
f{1}    = f1; f{2} = f2;

% use Valentin's function, at some point I should figure out the Matlab anova syntax
stats = rm_anova(y1(:), s(:), f);

% ========================================================= %
% BAYESIAN REPEATED MEASURES ANOVA
% ========================================================= %

% do Bayesian ANOVA to get Bayes Factors
statdat         = table;
statdat.DV      = y1(:);
statdat.subjnr  = s(:);
statdat.prevPupilBins = f2(:);
statdat.session = f1(:);
writetable(statdat, sprintf('%s/Data/CSV/sessionANOVAdat.csv', mypath));
system('/Library/Frameworks/R.framework/Resources/bin/R < BayesFactorANOVA_sessions.R --no-save');
statres = readtable(sprintf('%s/Data/CSV/sessionANOVAresults.csv', mypath)); % fetch results

fprintf('Session F(%d,%d) = %.2f, p = %.3f, Bf10 = %.3f \n', stats.f1.df(1), stats.f1.df(2), stats.f1.fstats, stats.f1.pvalue, statres.session);
fprintf('Interaction F(%d,%d) = %.2f, p = %.3f, Bf10 = %.3f  \n', stats.f1xf2.df(1), stats.f1xf2.df(2), stats.f1xf2.fstats, stats.f1xf2.pvalue, statres.pupil_session);

%% wite out for JASP format
statdat.prevPupilBins = categorical(statdat.prevPupilBins);
statdat.session = categorical(statdat.session);

test = unstack(statdat, 'DV', {'prevPupilBins'});
test = unstack(test, {'x1', 'x2', 'x3'}, 'session');
writetable(test, sprintf('%s/Data/CSV/JASP_AnovaCheck.csv', mypath));
