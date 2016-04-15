clear all; clc; close;

mypath = '/Users/anne/Data/pupilUncertainty';

nbins   = 3;
lag     = 1;
sj      = 6;
data    = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
data    = data(find(data.sessionnr > 1), :);
data.evidence = abs(data.motionstrength); % single-trial evidence strength is absolute motionenergy
data.resp(data.resp == -1) = 0; % logistic regression outcome
cnt = 1; % for subplots

flds = {'evidence',  'rt',};
flds = {'rt', 'evidence'};
% flds = {'evidence'};
for f = 1:length(flds),
    
    % pupil quantiles vs evidence
    subplot(4,3,cnt); hold on; cnt = cnt + 1;
    [binnedx, binnedy, stdx, stdy, rho, pval] = divideintobins( ...
        data.evidence(data.correct==1), data.(flds{f})(data.correct==1), 10);
    errorbar(binnedx, binnedy, stdy/10); xlabel('evidence'); ylabel(flds{f});
    axis tight; box off; axis square;
    
    resps = [0 1];
    for r = 1:2,
        % 3 bins of pupil dilation, response bias
        uncQs = quantile(data.(flds{f}), nbins-1);
        for u = 1:nbins,
            switch u
                case 1
                    trls = find(data.resp == resps(r) & data.correct == 1 & data.(flds{f}) < uncQs(u));
                case nbins
                    trls = find(data.resp == resps(r) & data.correct == 1 & data.(flds{f}) > uncQs(u-1));
                otherwise
                    trls = find(data.resp == resps(r) & data.correct == 1 & ...
                        data.(flds{f}) > uncQs(u-1) & data.(flds{f}) < uncQs(u));
            end
            assert(~isempty(trls));
            
            % with this selection, take the trials after that
            laggedtrls = trls+lag;
            % exclude trials at the end of the block
            if any(laggedtrls > size(data, 1)),
                trls(laggedtrls > size(data, 1)) = [];
                laggedtrls(laggedtrls > size(data, 1)) = [];
            end
            
            % remove trials that dont match in block nr
            removeTrls = data.blocknr(laggedtrls) ~= data.blocknr(trls);
            laggedtrls(removeTrls) = [];
            trls(removeTrls) = [];
            
            % fit logistic regression
            thisdat = data(laggedtrls, :);
            [b, dev, stats] = glmfit(thisdat.motionstrength, thisdat.resp, ...
                'binomial','link','logit');
            
            grandavg.(flds{f})(r, u) = exp(b(1)) ./ (1 + exp(b(1)));
            grandavg.RTcheck(r,u)    = median(thisdat.evidence);
        end
    end
    % collapse the 2 response identities
    resp1 = -1 * (squeeze(grandavg.(flds{f})(1,:)) - 0.5) + 0.5;
    resp2 = grandavg.(flds{f})(2,:);
    resp1 = grandavg.(flds{f})(1,:);
    
    grandavg.rep = resp1; %(resp1 + resp2) ./ 2;
    subplot(4,3,cnt); hold on; cnt = cnt + 1;
    
    grandavg.RT = mean(grandavg.RTcheck);
    plot(grandavg.RT); hold on; ylabel('evidence(t+1)'); xlabel(flds{f});
    axis square;
    
    subplot(4,3,cnt); hold on; cnt = cnt + 1;
    plot(resp1); hold on; plot(resp2); ylabel('rep'); xlabel(flds{f});
    axis square;
    % ylim([0.4 0.6]);
    title('correct trials');
    
end


