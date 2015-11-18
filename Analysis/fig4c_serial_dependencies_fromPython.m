%% this script requires that the intertrial toolbox in Python has been run for all participants
% read in the python-generated mat files and do some plots!

clear; close all; clc;
subjects = 26:27;
nlags = 7;
lags = 1:7;
colors = linspecer(8);

for sj = subjects,
    clf;
    
    % ============================================ %
    % ======= model without pupil term =========== %
    % ============================================ %
    
    load(sprintf('~/Dropbox/code/serial-dependencies/sim_backup/2ifc_sj%02d.txtresults.mat', sj));
    load(sprintf('~/Dropbox/code/serial-dependencies/sim_backup/2ifc_sj%02d.txtdata.mat', sj));
    
    % get the fitted weights
    hf0 = model_w_hist.hf0+1;
    hlen = size(h, 2);
    
    respw = model_w_hist.w(hf0:hf0+hlen-1);
    stimw = model_w_hist.w(hf0+hlen:hf0+hlen*2-1);
    
    % project back into lag space
    stimw = h * stimw';
    respw = h * respw';
    
    % multiply by the mean slope to get more meaningful units
    slopeidx = 2:hf0-1;
    % alpha = mean(model_w_hist.w(slopeidx));
    %stimw = stimw ./ alpha;
    %respw = respw ./ alpha;
    
    % for bootstrapped values, also multiply with the bootstrapped slope
    bootstrap_corr(:,1:nlags*2)  = bsxfun(@times, bootstrap(:, 1:nlags*2), bootstrap(:, end-1));
    
    % also get error bars, multiply bootstrapped values with slope too
    respboot = bootstrap_corr(:, 1:nlags);
    stimboot = bootstrap_corr(:, nlags+1:nlags*2);
    
    respwci = prctile(respboot, [2.5, 97.5])';
    respwci(:, 1) = respw - respwci(:, 1); respwci(:, 2) = respwci(:, 2) - respw;
    stimwci = prctile(stimboot, [2.5, 97.5])';
    stimwci(:, 1) = stimw - stimwci(:, 1); stimwci(:, 2) = stimwci(:, 2) - stimw;
    
    % ============================================ %
    % 1. plot the history kernels for resp and stim (blue/yellow as Fr?nd)
    
    subplot(3,3,1); hold on;
    bh = boundedline(lags, stimw, stimwci, lags, respw, respwci, 'cmap', colors([2 4], :), 'alpha');
    legend(bh, 'stimulus', 'response'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    ylabel('No pupil');
    
    % ============================================ %
    % 2. same thing, but for correct and error (green/red)
    correctw = stimw + respw;
    errorw   = -stimw + respw;
    correctboot = stimboot + respboot;
    errorboot = -stimboot + respboot;
    
    correctwci = prctile(correctboot, [2.5, 97.5])';
    correctwci(:, 1) = correctw - correctwci(:, 1); correctwci(:, 2) = correctwci(:, 2) - correctw;
    errorwci = prctile(errorboot, [2.5, 97.5])';
    errorwci(:, 1) = errorw - errorwci(:, 1); errorwci(:, 2) = errorwci(:, 2) - errorw;
    
    subplot(3,3,2);
    bh =  boundedline(lags, correctw, correctwci, lags, errorw, errorwci, 'cmap', colors([3 1], :), 'alpha');
    legend(bh, 'correct', 'error'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    
    % ============================================ %
    % 3. plot the psychometric functions for each session
    % get the data structure for this
    
    clear yfit
    for cond = 1:6,
        bias = model_w_hist.w(1);
        slope = model_w_hist.w(1+cond);
        
        % glmval and plot the curve
        x = [-.3:0.001:0.3];
        yfit(:, cond) = glmval([bias slope]', x,'logit');
    end
    subplot(3,3,3);
    set(groot,'defaultAxesColorOrder',parula(10))
    plot(x, yfit);
    box off; xlim([-0.3 0.3]);
    xlabel('Stimulus'); ylabel('P(respA) over sessions');
    
    % ============================================ %
    % ======= model WITH pupil term =========== %
    % ============================================ %
    
    load(sprintf('~/Dropbox/code/serial-dependencies/sim_backup/2ifc_pupil_sj%02d.txtresults.mat', sj));
    load(sprintf('~/Dropbox/code/serial-dependencies/sim_backup/2ifc_pupil_sj%02d.txtdata.mat', sj));
    
    % get the fitted weights
    hf0 = model_w_hist.hf0+1;
    hlen = size(h, 2);
    
    respw = model_w_hist.w(hf0:hf0+hlen-1);
    stimw = model_w_hist.w(hf0+hlen:hf0+hlen*2-1);
    
    % project back into lag space
    stimw = h * stimw';
    respw = h * respw';
    
    % for bootstrapped values, also multiply with the bootstrapped slope
    bootstrap_corr(:,1:size(bootstrap,2)-2)  = bsxfun(@times, bootstrap(:, 1:end-2), bootstrap(:, end-1));
    
    % also get error bars, multiply bootstrapped values with slope too
    respboot = bootstrap_corr(:, 1:nlags);
    stimboot = bootstrap_corr(:, nlags+1:nlags*2);
    
    respwci = prctile(respboot, [2.5, 97.5])';
    respwci(:, 1) = respw - respwci(:, 1); respwci(:, 2) = respwci(:, 2) - respw;
    stimwci = prctile(stimboot, [2.5, 97.5])';
    stimwci(:, 1) = stimw - stimwci(:, 1); stimwci(:, 2) = stimwci(:, 2) - stimw;
    
    % ============================================ %
    % 1. plot the history kernels for resp and stim (blue/yellow as Fr?nd)
    subplot(3,3,4); hold on;
    bh = boundedline(lags, stimw, stimwci, lags, respw, respwci, 'cmap', colors([2 4], :), 'alpha');
    legend(bh, 'stimulus', 'response'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    ylabel('With pupil');
    
    % ============================================ %
    % 2. same thing, but for correct and error (green/red)
    correctw = stimw + respw;
    errorw   = -stimw + respw;
    correctboot = stimboot + respboot;
    errorboot = -stimboot + respboot;
    
    correctwci = prctile(correctboot, [2.5, 97.5])';
    correctwci(:, 1) = correctw - correctwci(:, 1); correctwci(:, 2) = correctwci(:, 2) - correctw;
    errorwci = prctile(errorboot, [2.5, 97.5])';
    errorwci(:, 1) = errorw - errorwci(:, 1); errorwci(:, 2) = errorwci(:, 2) - errorw;
    
    subplot(3,3,5);
    bh =  boundedline(lags, correctw, correctwci, lags, errorw, errorwci, 'cmap', colors([3 1], :), 'alpha');
    legend(bh, 'correct', 'error'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    
    % ============================================ %
    % 3. the pure pupil weights
    
    pupilw = model_w_hist.w(hf0+hlen*2:hf0+hlen*3-1);
    pupilw = h * pupilw';
    pupilboot = bootstrap_corr(:, nlags*2+1:nlags*3);
    
    pupilwci = prctile(pupilboot, [2.5, 97.5])';
    pupilwci(:, 1) = pupilw - pupilwci(:, 1); pupilwci(:, 2) = pupilwci(:, 2) - pupilw;
    
    subplot(3,3,6); hold on;
    bh = boundedline(lags, pupilw, pupilwci, 'cmap', colors(5, :));
    legend(bh, 'pupil'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    
    % ============================================ %
    % 4. pupil * history weights
    
    pupilrespw = model_w_hist.w(hf0+hlen*3:hf0+hlen*4-1);
    pupilstimw = model_w_hist.w(hf0+hlen*4:hf0+hlen*5-1);
    
    % project back into lag space
    pupilrespw = h * pupilrespw';
    pupilstimw = h * pupilstimw';

    % also get error bars, multiply bootstrapped values with slope too
    pupilrespboot = bootstrap_corr(:, nlags*3+1:nlags*4);
    pupilstimboot = bootstrap_corr(:, nlags*4+1:nlags*5);
    
    pupilrespwci = prctile(pupilrespboot, [2.5, 97.5])';
    pupilrespwci(:, 1) = pupilrespw - pupilrespwci(:, 1); pupilrespwci(:, 2) = pupilrespwci(:, 2) - pupilrespw;
    pupilstimwci = prctile(pupilstimboot, [2.5, 97.5])';
    pupilstimwci(:, 1) = pupilstimw - pupilstimwci(:, 1); pupilstimwci(:, 2) = pupilstimwci(:, 2) - pupilstimw;
    
    % 1. plot the history kernels for resp and stim (blue/yellow as Fr?nd)
    subplot(3,3,7); hold on;
    bh = boundedline(lags, pupilstimw, pupilstimwci, lags, pupilrespw, pupilrespwci, 'cmap', colors([2 4], :), 'alpha');
    legend(bh, 'pupil*stimulus', 'pupil*response'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    ylabel('Interaction');
    
    % 2. same thing, but for correct and error (green/red)
    pupilcorrectw = pupilstimw + pupilrespw;
    pupilerrorw   = -pupilstimw + pupilrespw;
    pupilcorrectboot = pupilstimboot + pupilrespboot;
    pupilerrorboot = -pupilstimboot + pupilrespboot;
    
    pupilcorrectwci = prctile(pupilcorrectboot, [2.5, 97.5])';
    pupilcorrectwci(:, 1) = pupilcorrectw - pupilcorrectwci(:, 1); pupilcorrectwci(:, 2) = pupilcorrectwci(:, 2) - pupilcorrectw;
    pupilerrorwci = prctile(pupilerrorboot, [2.5, 97.5])';
    pupilerrorwci(:, 1) = pupilerrorw - pupilerrorwci(:, 1); pupilerrorwci(:, 2) = pupilerrorwci(:, 2) - pupilerrorw;
    
    subplot(3,3,8);
    bh =  boundedline(lags, pupilcorrectw, pupilcorrectwci, lags, pupilerrorw, pupilerrorwci, 'cmap', colors([3 1], :), 'alpha');
    legend(bh, 'pupil*correct', 'pupil*error'); legend boxoff;
    xlim([0.5 nlags+0.5]);
    
    % ============================================ %
    % ======= model WITH pupil term =========== %
    % ============================================ %
    
    print(gcf, '-dpdf', sprintf('~/Dropbox/Figures/sequential/sequential_withpupil_sj%02d.pdf', sj));

end

