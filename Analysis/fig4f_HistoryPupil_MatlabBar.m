function fig4f_HistoryPupil_MatlabBar(lagGroups, whichmodulator)

if ~exist('lagGroups', 'var'), lagGroups = 1; end

% ============================================ %
% dont use Fruend but my own regressors
% ============================================ %

subjects = 1:27;
grandavg.betas = nan(27, 3, 7);
grandavg.pval = nan(27, 3, 7);

for sj = unique(subjects),
    data = readtable(sprintf('~/Data/pupilUncertainty/CSV/2ifc_data_sj%02d.csv', sj));
    
    % get the variables we need
    resp = data.resp; resp(resp == -1) = 0; % predict response identity
    motionstrength = zscore(data.motionstrength);
    prevPupil      = circshift(zscore(data.decision_pupil), 1);
    prevResp       = circshift(data.resp, 1);
    prevStim       = circshift(data.stim, 1);
    
    % make design matrix - intercept will be added automatically
    designM = [motionstrength prevResp prevResp.*prevPupil];
    xlabs = {'bias', 'stim', 'prevResp','prevResp*prevPupil'};
    
    % don't use trials that are at the beginning of each block
    trlDif = [0; diff(data.trialnr)];
    
    removeTrls = false(size(trlDif));
    removeTrls(trlDif < 1) = true;
    removeTrls(trlDif > 1) = true;
    removeTrls(find(trlDif > 1) + 1) = true;
    
    % in these trials, the history wont be able to predict the response
    designM(removeTrls==1, 2:end) = 0;
     
    % do this separately for error, correct and all trials together
    cors = [0 1];
    for c = 1:3,
        correctPrev = circshift(data.correct, 1);
        if c == 3,
            trls = 1:length(resp);
        else
            trls = find(correctPrev == cors(c));
        end

        mdl  = fitglm(designM(trls, :), resp(trls), ...
            'distr', 'binomial', 'link', 'logit');
        grandavg.betas(sj, c, 1:length(mdl.Coefficients.Estimate)) = mdl.Coefficients.Estimate;
        grandavg.pval(sj, c, 1:length(mdl.Coefficients.pValue))  = mdl.Coefficients.pValue;
        
        grandavg.dev(sj, 3) = mdl.Deviance;
        
        if c == 3, % only on all trials
            
            % get the deviance for the null model without history
            mdl  = fitglm(designM(trls, 1), resp(trls), ...
                'distr', 'binomial', 'link', 'logit');
            grandavg.dev(sj, 1) = mdl.Deviance;
            
            % model with  history
            mdl  = fitglm(designM(trls, 1:2), resp(trls), ...
                'distr', 'binomial', 'link', 'logit');
            grandavg.dev(sj, 2) = mdl.Deviance;
            
        end
    end
end

% ============================================ %
% combine lags
% ============================================ %

groupedDat.response = grandavg.betas(:, 3, 3);
groupedDat.response_pupil = grandavg.betas(:, 3, 4);
groupedDat.correct = grandavg.betas(:, 2, 3);
groupedDat.correct_pupil = grandavg.betas(:, 2, 4);
groupedDat.incorrect = grandavg.betas(:, 1, 3);
groupedDat.incorrect_pupil = grandavg.betas(:, 1, 4);

% ============================================ %
% barweb matrix
% ============================================ %

lagGroups = 1;
posRespSj = find(mean(groupedDat.response(:, lagGroups), 2) > 0);
negRespSj = find(mean(groupedDat.response(:, lagGroups), 2) < 0);

load(sprintf('~/Data/pupilUncertainty/GrandAverage/historyweights_%s.mat', 'pupil'));
%posRespSj = find(mean(dat.response(:, lagGroups), 2) > 0);
%negRespSj = find(mean(dat.response(:, lagGroups), 2) < 0);

colors = linspecer(9);

bwMat = cat(3, [groupedDat.response groupedDat.correct groupedDat.incorrect], ...
    [groupedDat.response_pupil groupedDat.correct_pupil groupedDat.incorrect_pupil]);

% sigstars
pvals.posRespSj = squeeze(cat(3, [randtest1d(groupedDat.response(posRespSj)) randtest1d(groupedDat.correct(posRespSj)) randtest1d(groupedDat.incorrect(posRespSj))], ...
    [randtest1d(groupedDat.response_pupil(posRespSj)) randtest1d(groupedDat.correct_pupil(posRespSj)) randtest1d(groupedDat.incorrect_pupil(posRespSj))]))';
pvals.negRespSj = squeeze(cat(3, [randtest1d(groupedDat.response(negRespSj)) randtest1d(groupedDat.correct(negRespSj)) randtest1d(groupedDat.incorrect(negRespSj))], ...
    [randtest1d(groupedDat.response_pupil(negRespSj)) randtest1d(groupedDat.correct_pupil(negRespSj)) randtest1d(groupedDat.incorrect_pupil(negRespSj))]))';

subplot(4,4,10); hold on;
barwebQH(squeeze(nanmean(bwMat(posRespSj, :, :)))', squeeze(nanstd(bwMat(posRespSj, :, :)))' ./sqrt(length(posRespSj)), ...
    [], 0.5, colors([2 3 1], :), pvals.posRespSj);
set(gca, 'xtick', 1:2, 'xticklabel', {'Resp-1', 'Resp-1*Pup-1'});
%axis tight; 
ylim([-.15 0.3]);
axis square; xlim([0.5 2.5]);
title('Repeaters', 'color', colors(8,:));

subplot(4,4,11); hold on;
barwebQH(squeeze(nanmean(bwMat(negRespSj, :, :)))', squeeze(nanstd(bwMat(negRespSj, :, :)))' ./sqrt(length(negRespSj)), ...
    [], 0.5, colors([2 3 1], :), pvals.negRespSj);
title('Switchers', 'color', colors(9,:));
axis square; xlim([0.5 2.5]);
ylim([-.35 0.1]);
set(gca, 'xtick', 1:2, 'xticklabel',{'Resp-1', 'Resp-1*Pup-1'});

if ~exist('whichmodulator', 'var'), whichmodulator = 'pupil'; end

print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/uncertainty/fig4e_historyPupil_%s.pdf', whichmodulator));
