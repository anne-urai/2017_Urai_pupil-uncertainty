function grandavg = postPupilBehaviour(whichmodulator, nbins, correctness)
global mypath;

if ~exist('whichmodulator', 'var'); whichmodulator = 'pupil'; end
if ~exist('correctness', 'var'); correctness = []; end
if ~exist('nbins', 'var'); nbins = 3; end
lag = 1; % look at 1 trial in the past

warning('error', 'stats:glmfit:PerfectSeparation');

% =========================================== %

subjects = 1:27;
clear grandavg;

for sj = unique(subjects),
    
    % get data
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    data = data((data.sessionnr > 1), :); % get rid of residual learning effects
    
    % determine what we will bin by
    switch whichmodulator
        case 'fb+decpupil'
            whichMod = 'feedback_pupil';
            data.feedback_pupil = projectout(data.feedback_pupil, data.decision_pupil);
        case 'dec+fbpupil'                    
            whichMod = 'decision_pupil';
            data.decision_pupil = projectout(data.decision_pupil, data.feedback_pupil);
        case 'pupil',
            whichMod = 'decision_pupil';
            data.decision_pupil = projectout(data.decision_pupil, data.rtNorm); % take out RT
        case 'rt'
            whichMod     = 'rtNorm';
            data.rtNorm = projectout(data.rtNorm, data.decision_pupil);
        case 'evidence'
            % single-trial evidence strength is absolute motionenergy
            data.evidence = abs(data.motionstrength);
            % make sure to use the same range of evidence for each session
            for session = unique(data.sessionnr)',
                data.evidence(data.sessionnr == session) = ...
                    zscore(data.evidence(data.sessionnr == session));
            end
        otherwise
            whichMod = whichmodulator;
    end
    
    % outcome vector need to be 0 1 for logistic regression
    data.resp(data.resp == -1) = 0;
    
    % get an overall logistic fit for normalization
    grandavg.overallLogistic(sj, :) = glmfit(data.motionstrength, ...
        data.resp, 'binomial','link','logit');
    
    % split into quantiles of the modulator
    if isempty(correctness),
        uncQs = quantile(data.(whichMod), nbins - 1);
    else
        uncQs = quantile(data.(whichMod)(data.correct == correctness), nbins - 1);
    end
    
    % =========================================== %
    % uncertainty bins
    % =========================================== %

    for u = 1:nbins,
        
        % find the trials in this bin
        if isempty(correctness), % take all trials
            switch u
                case 1
                    trls = find(data.(whichMod) <= uncQs(u));
                case nbins % last one
                    trls = find(data.(whichMod) > uncQs(u-1));
                otherwise
                    trls = find(data.(whichMod) > uncQs(u-1) & data.(whichMod) <= uncQs(u));
            end
            
        else % either correct or error trials
            switch u
                case 1
                    trls = find(data.correct == correctness & data.(whichMod) <= uncQs(u));
                case nbins % last one
                    trls = find(data.correct == correctness & data.(whichMod) > uncQs(u-1));
                otherwise
                    trls = find(data.correct == correctness & ...
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
        laggedtrls(data.blocknr(laggedtrls) ~= data.blocknr(trls)) = [];
        
        % use these data
        thisdat = data(laggedtrls, :);
        
        % =========================================== %
        % overall logistic
        % =========================================== %
        
        try
            b = glmfit(thisdat.motionstrength, thisdat.resp, ...
                'binomial','link','logit');
        catch % if there is perfect separation, this will catch (see error warning above)
            warning('putting nan in betas!');
            b = nan(size(b));
        end
        
        grandavg.logistic(sj, u, :) = b;
        
        % =========================================== %
        % overall RT, median
        % =========================================== %
        
        grandavg.RT(sj, u)       = nanmedian(thisdat.rt);
        grandavg.accuracy(sj, u) = nanmean(thisdat.correct);
        
        % =========================================== %
        % history-dependent logistic
        % =========================================== %

        % previous response
        resps = [0 1];
        for r = 1:2,
            
            % take a further subset of trials with previous response 1 or 0
            laggedresptrls  = find(data.resp == resps(r));
            laggedresptrls  = laggedresptrls + 1; % go one trial later
            laggedresptrls  = intersect(laggedtrls, laggedresptrls);
            thisdat         = data(laggedresptrls, :);
            
            try
                b = glmfit(thisdat.motionstrength, thisdat.resp, ...
                    'binomial','link','logit');
            catch % if there is perfect separation, this will catch (see error warning above)
                warning('putting nan in betas!');
                b = nan(size(b));
            end
            
            % save betas
            grandavg.logisticHistory(sj, r, u, :) = b;
            
            % =========================================== %
            % Post-error slowing, Dutilh et al. 2012
            % =========================================== %
            
            errortrls = intersect(trls, find(data.correct == 0));
            
            % remove trls that are not in a continuous sequence
            errortrls(errortrls > size(data, 1)) = [];
            errortrls(errortrls < 2) = [];
            errortrls((data.trialnr(errortrls) - data.trialnr(errortrls - 1)) ~= 1) = [];
            errortrls((data.trialnr(errortrls) - data.trialnr(errortrls + 1)) ~= -1) = [];
            
            % use only those trials where both the pre-error and the
            % post-error are correct
            errortrls(data.correct(errortrls-1) == 0) = [];
            errortrls(data.correct(errortrls+1) == 0) = [];
            
            % see very clear image http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4807057/figure/pone.0151763.g001/
            grandavg.pes(sj, u) =  mean(data.rt(errortrls + 1) - data.rt(errortrls - 1));
            
        end % previous resp
    end % uncertainty bin
end % sj

% ========================================================= %
% transform history-dependent bias into repetition bias
% ========================================================= %

resp1 = -1 * squeeze(grandavg.logisticHistory(:, 1, :, 1)); % take the 1st beta weigth only, ignore history-dependent slope
resp2 =      squeeze(grandavg.logisticHistory(:, 2, :, 1));

% since this is centred at 0.5, treat it that way
grandavg.repetition = (resp1 + resp2) ./ 2;

% ========================================================= %
% transform from log(odds) to probability
% ========================================================= %

logOdds2Probability = @(x) exp(x) ./ (1 + exp(x));

grandavg.repetition = logOdds2Probability(grandavg.repetition);
% also for saying 1 over 2 - get rid of the slope
grandavg.logisticHistory = logOdds2Probability(grandavg.logisticHistory(:, :, :, 1));

% ========================================================= %
% pull out sensitivity and absolute side bias
% ========================================================= %

grandavg.sensitivity    = grandavg.logistic(:, :, 2);
grandavg.signedBias     = grandavg.logistic(:, :, 1);
grandavg.absoluteBias   = abs(grandavg.logistic(:, :, 1));

end
