function a6_retrieveDataFromPython(whichmodulator)
%% this script requires that the intertrial toolbox in Python has been run for all participants
% read in the python-generated mat files and do some plots!

global mypath;

clear; close all; clc;
subjects = 1:27;
nlags = 7;
lags = 1:7;
colors = linspecer(8);

% preallocate
dat.response = nan(27, nlags);
dat.responseCI = nan(27, nlags, 2);
dat.stimulus = nan(27, nlags);
dat.stimulusCI = nan(27, nlags, 2);
dat.correct = nan(27, nlags);
dat.incorrect = nan(27, nlags);
dat.pupil = nan(27, nlags);
dat.response_pupil = nan(27, nlags);
dat.stimulus_pupil = nan(27, nlags);
dat.correct_pupil = nan(27, nlags);
dat.incorrect_pupil = nan(27, nlags);
dat.historyPval = nan(27, 1);
dat.modulationPval = nan(27, 1);

for sj = subjects,
    clf;
    
    % ============================================ %
    % ======= model WITH pupil term =========== %
    % ============================================ %
    try
        load(sprintf('%s/Data/serialmodel/2ifc_%s_sj%02d.txtresults.mat', mypath, whichmodulator, sj));
        load(sprintf('%s/Data/serialmodel/2ifc_%s_sj%02d.txtdata.mat', mypath, whichmodulator, sj));
    catch
        warning('skipping participant %02d', sj);
        continue;
    end
    
    switch whichmodulator
        case 'plain'
            model_h_mod = model_w_hist;
    end
    
    % get the fitted weights
    hf0 = model_h_mod.hf0+1;
    model_h_mod.w = model_h_mod.w(:);
    hlen = size(h, 2);
    
    respw = model_h_mod.w(hf0:hf0+hlen-1);
    stimw = model_h_mod.w(hf0+hlen:hf0+hlen*2-1);
    
    % project back into lag space
    stimw = h * stimw;
    respw = h * respw;
    
    % for bootstrapped values, also multiply with the bootstrapped slope
    bootstrap_corr(:,1:size(bootstrap,2)-2)  = bsxfun(@times, bootstrap(:, 1:end-2), bootstrap(:, end-1));
    
    % also get error bars, multiply bootstrapped values with slope too
    respboot    = bootstrap_corr(:, 1:nlags);
    stimboot    = bootstrap_corr(:, nlags+1:nlags*2);
    alpha       = 1 - 0.68; % should cover 1 std of the distribution
    respwci     = prctile(respboot, [100*alpha/2,100*(1-alpha/2)])';
    dat.responseCI(sj, :, :) = respwci;
    
    respwci(:, 1) = respw - respwci(:, 1); respwci(:, 2) = respwci(:, 2) - respw; % relative error
    stimwci       = prctile(stimboot, [100*alpha/2,100*(1-alpha/2)])';
    dat.stimulusCI(sj, :, :) = stimwci;
    
    stimwci(:, 1) = stimw - stimwci(:, 1); stimwci(:, 2) = stimwci(:, 2) - stimw;
    
    % ============================================ %
    % 1. plot the history kernels for resp and stim (blue/yellow as Fr?nd)
    subplot(3,3,1); hold on;
    bh = boundedline(lags, stimw, stimwci, lags, respw, respwci, 'cmap', colors([2 4], :), 'alpha');
    legend(bh, 'stimulus', 'response'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    ylabel('History kernels');
    
    dat.response(sj, :) = respw;
    dat.stimulus(sj, :) = stimw;
    
    % ============================================ %
    % 2. same thing, but for correct and error (green/red)
    correctw = stimw + respw;
    errorw   = -stimw + respw;
    correctboot = stimboot + respboot;
    errorboot = -stimboot + respboot;
    
    % correct - error = stim + resp --stim +- resp = 2stim
    % error - correct = -stim + resp -resp - stim = -2stim
    
    correctwci = prctile(correctboot, [2.5, 97.5])';
    correctwci(:, 1) = correctw - correctwci(:, 1); correctwci(:, 2) = correctwci(:, 2) - correctw;
    errorwci = prctile(errorboot, [2.5, 97.5])';
    errorwci(:, 1) = errorw - errorwci(:, 1); errorwci(:, 2) = errorwci(:, 2) - errorw;
    
    subplot(3,3,2);
    bh =  boundedline(lags, correctw, correctwci, lags, errorw, errorwci, 'cmap', colors([3 1], :), 'alpha');
    legend(bh, 'correct', 'error'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    
    dat.correct(sj, :) = correctw;
    dat.incorrect(sj, :) = errorw;
    
    % ============================================ %
    % 3. plot the psychometric functions for each session
    % get the data structure for this
    % ============================================ %
    
    clear yfit
    for cond = 1:6,
        bias = model_h_mod.w(1);
        slope = model_h_mod.w(1+cond);
        
        % glmval and plot the curve
        x = [-.3:0.001:0.3];
        yfit(:, cond) = glmval([bias slope]', x,'logit');
    end
    subplot(3,3,3);
    set(groot,'defaultAxesColorOrder',parula(10))
    plot(x, yfit);
    box off; xlim([-0.3 0.3]);
    xlabel('Stimulus'); ylabel('P(respA) over sessions');
    
    switch whichmodulator
        case 'plain'
        otherwise % only do when we ran a modulatory model
            
            % ============================================ %
            % 3. the pure pupil weights
            % ============================================ %
            
            pupilw = model_h_mod.w(hf0+hlen*2:hf0+hlen*3-1);
            pupilw = h * pupilw;
            pupilboot = bootstrap_corr(:, nlags*2+1:nlags*3);
            
            pupilwci = prctile(pupilboot, [2.5, 97.5])';
            pupilwci(:, 1) = pupilw - pupilwci(:, 1); pupilwci(:, 2) = pupilwci(:, 2) - pupilw;
            
            subplot(3,3,6); hold on;
            bh = boundedline(lags, pupilw, pupilwci, 'cmap', colors(5, :));
            legend(bh, whichmodulator); legend boxoff;
            xlim([0.5 nlags+0.5]);
            
            dat.pupil(sj, :) = pupilw;
            
            % ============================================ %
            % 4. pupil * history weights
            % ============================================ %
            
            pupilrespw = model_h_mod.w(hf0+hlen*3:hf0+hlen*4-1);
            pupilstimw = model_h_mod.w(hf0+hlen*4:hf0+hlen*5-1);
            
            % project back into lag space
            pupilrespw = h * pupilrespw;
            pupilstimw = h * pupilstimw;
            
            % also get error bars, multiply bootstrapped values with slope too
            pupilrespboot = bootstrap_corr(:, nlags*3+1:nlags*4);
            pupilstimboot = bootstrap_corr(:, nlags*4+1:nlags*5);
            
            pupilrespwci = prctile(pupilrespboot, [2.5, 97.5])';
            pupilrespwci(:, 1) = pupilrespw - pupilrespwci(:, 1); pupilrespwci(:, 2) = pupilrespwci(:, 2) - pupilrespw;
            pupilstimwci = prctile(pupilstimboot, [2.5, 97.5])';
            pupilstimwci(:, 1) = pupilstimw - pupilstimwci(:, 1); pupilstimwci(:, 2) = pupilstimwci(:, 2) - pupilstimw;
            
            % 1. plot the history kernels for resp and stim (blue/yellow as Fr?nd)
            subplot(3,3,4); hold on;
            bh = boundedline(lags, pupilstimw, pupilstimwci, lags, pupilrespw, pupilrespwci, 'cmap', colors([2 4], :), 'alpha');
            legend(bh, [whichmodulator '*stimulus'], [whichmodulator '*response']); legend boxoff;
            xlim([0.5 nlags+0.5]);
            ylabel('Interaction');
            
            dat.response_pupil(sj, :) = pupilrespw;
            dat.stimulus_pupil(sj, :) = pupilstimw;
            
            % 2. same thing, but for correct and error (green/red)
            pupilcorrectw = pupilstimw + pupilrespw;
            pupilerrorw   = -pupilstimw + pupilrespw;
            pupilcorrectboot = pupilstimboot + pupilrespboot;
            pupilerrorboot = -pupilstimboot + pupilrespboot;
            
            pupilcorrectwci = prctile(pupilcorrectboot, [2.5, 97.5])';
            pupilcorrectwci(:, 1) = pupilcorrectw - pupilcorrectwci(:, 1); pupilcorrectwci(:, 2) = pupilcorrectwci(:, 2) - pupilcorrectw;
            pupilerrorwci = prctile(pupilerrorboot, [2.5, 97.5])';
            pupilerrorwci(:, 1) = pupilerrorw - pupilerrorwci(:, 1); pupilerrorwci(:, 2) = pupilerrorwci(:, 2) - pupilerrorw;
            
            subplot(3,3,5);
            bh =  boundedline(lags, pupilcorrectw, pupilcorrectwci, lags, pupilerrorw, pupilerrorwci, 'cmap', colors([3 1], :), 'alpha');
            legend(bh, [whichmodulator '*correct'], [whichmodulator '*error']); legend boxoff;
            xlim([0.5 nlags+0.5]);
            
            dat.correct_pupil(sj, :) = pupilcorrectw;
            dat.incorrect_pupil(sj, :) = pupilerrorw;
            
            suplabel(sprintf('P%02d', sj), 't');
            
            % ============================================ %
            % 5. permutation testing
            % ============================================ %
            
            subplot(3,3,7);
            histogram(permutation_wh(:, 1), 'EdgeColor', 'none', 'facecolor', [0.8 0.8 0.8]); hold on;
            axis tight; box off;
            pval = length(find(permutation_wh(:, 1) > model_w_hist.loglikelihood)) ./ length(permutation_wh(:, 1));
            prc95 = prctile(permutation_wh(:, 1), 97.5);
            plot([prc95 prc95], [0 max(get(gca, 'ylim'))], 'k');
            plot([model_w_hist.loglikelihood model_w_hist.loglikelihood], [0 max(get(gca, 'ylim'))], 'r');
            xlabel('logLikelihood history');
            dat.historyPval(sj) = pval;
            
            % now for the pupil
            subplot(3,3,8);
            histogram(permutation_hmod(:, 1), 'EdgeColor', 'none', 'facecolor', [0.8 0.8 0.8]); hold on;
            box off;
            pval = length(find(permutation_hmod(:, 1) > model_h_mod.loglikelihood)) ./ length(permutation_hmod(:, 1));
            if pval < 0.05, color = 'r'; else color = 'b'; end
            prc95 = prctile(permutation_hmod(:, 1), 97.5);
            pl(1) = plot([prc95 prc95], [0 max(get(gca, 'ylim'))], 'k');
            pl(2) = plot([model_w_hist.loglikelihood model_w_hist.loglikelihood], [0 max(get(gca, 'ylim'))], 'r'); % history loglik
            pl(3) = plot([model_h_mod.loglikelihood model_h_mod.loglikelihood], [0 max(get(gca, 'ylim'))], 'b');
            xlabel('logLikelihood modulation');
            dat.modulationPval(sj) = pval;
            
            lh = legend(pl, {'95th percentile', 'loglik history-only', 'loglik history*pupil'});
            sp = subplot(3,3,9);
            spos = get(sp, 'position'); set(lh, 'position', spos, 'box', 'off'); axis off
    end
    print(gcf, '-dpdf', sprintf('%s/Figures/serial/historykernels_%s_sj%02d.pdf', mypath, whichmodulator, sj));
    
end

% save for group plots
savefast(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator), 'dat');
close all;
end