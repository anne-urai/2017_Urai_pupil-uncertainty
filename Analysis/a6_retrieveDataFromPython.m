function a6_retrieveDataFromPython(whichmodulator)
%% this script requires that the intertrial toolbox in Python has been run for all participants
% read in the python-generated mat files and do some plots!

global mypath;
subjects = 1:27;
nlags = 7;
disp(whichmodulator);

for sj = subjects,
    % ============================================ %
    % ======= model WITH pupil term =========== %
    % ============================================ %
    
    load(sprintf('%s/Data/serialmodel/2ifc_%s_sj%02d.txtresults.mat', mypath, whichmodulator, sj));
    load(sprintf('%s/Data/serialmodel/2ifc_%s_sj%02d.txtdata.mat', mypath, whichmodulator, sj));
    
    switch whichmodulator
        case {'plain', 'plainCoh'} % no modulator
            model_h_mod = model_w_hist;
            weights1 = {'response', 'stimulus'};
            weights2 = {''};
        case {'plain_session2', 'plain_session3', 'plain_session4', ...
                'plain_session5', 'plain_session6'}; % by session
            model_h_mod = model_w_hist;
            weights1 = {'response', 'stimulus'};
            weights2 = {''};
        case {'fbpupil', 'pupil', 'rt'}, % single modulator
            weights1 = {'response', 'stimulus', ...
                'pupil', 'response_pupil', 'stimulus_pupil'};
            weights2 = {'', '_pupil'};
        otherwise % two modulators
            weights1 = {'response', 'stimulus', ...
                'pupil', 'response_pupil', 'stimulus_pupil', ...
                'rt', 'response_rt', 'stimulus_rt'};
            weights2 = {'', '_pupil', '_rt'};
    end
    
    % get the fitted weights
    hf0 = model_h_mod.hf0+1;
    model_h_mod.w = model_h_mod.w(:);
    hlen = size(h, 2);
    
    % for bootstrapped values, also multiply with the bootstrapped slope
    bootstrap_corr(:,1:size(bootstrap,2)-2)  = bsxfun(@times, bootstrap(:, 1:end-2), bootstrap(:, end-1));

    cnt = 0;
    for w = 1:length(weights1),
        
        % project back into lag space
        thisw       = model_h_mod.w(hf0+hlen*cnt:hf0+hlen*(cnt+1)-1);
        thisw       = h * thisw;
        dat.(weights1{w})(sj, :) = thisw;
        
        % get errorbars
        thisboot    = bootstrap_corr(:, nlags*cnt+1:nlags*(cnt+1));
        alpha       = 1 - 0.68; % should cover 1 std of the distribution
        thiswci     = prctile(thisboot, [100*alpha/2,100*(1-alpha/2)])';
        dat.([weights1{w} 'CI'])(sj, :, :) = thiswci;
        cnt         = cnt  + 1;
    end
    
    % change into correct and error from stimulus and response
    % see Fruend et al, supplement A4
    % note: these are identical to reward and punishment weights 
    % from Busse et al, Abrahamyan et al

    cnt = 0;
    for w = 1:length(weights2),
        dat.(['correct' weights2{w}])(sj, :) = ...
            dat.(['stimulus' weights2{w}])(sj, :) + dat.(['response' weights2{w}])(sj, :);
        dat.(['incorrect' weights2{w}])(sj, :) = ...
            -dat.(['stimulus' weights2{w}])(sj, :) + dat.(['response' weights2{w}])(sj, :);
        
        % make correct and error CIs by recoding each bootstrap
        boot_choice    = bootstrap_corr(:, nlags*cnt+1:nlags*(cnt+1)); cnt = cnt + 1;
        boot_stim      = bootstrap_corr(:, nlags*cnt+1:nlags*(cnt+1)); cnt = cnt + 2;
        
        % combine them
        boot_correct   = boot_choice + boot_stim;
        boot_incorrect = boot_choice - boot_stim;

        % CIs in absolute bound values
        dat.(['correct' weights2{w} 'CI'])(sj, :, :) = prctile(boot_correct, [100*alpha/2,100*(1-alpha/2)])';
        dat.(['incorrect' weights2{w} 'CI'])(sj, :, :) = prctile(boot_incorrect, [100*alpha/2,100*(1-alpha/2)])';
        
    end

    % add the p value for permutation with history
    dat.pvalue(sj) = length(find(model_w_hist.loglikelihood < permutation_wh(:, 1))) ./ length(permutation_wh(:, 1));

    % if possible, get the variance explained
    try
        dat.variance.stimuli(sj, :)   = model_w_hist.stimuli;
        dat.variance.explained(sj, :) = model_w_hist.variance_explained;
    end
end

dat.pvalue = dat.pvalue';
switch whichmodulator
    case 'plain'
        % how many alternators have significant history couplings?
        altSignfic = length(find(dat.pvalue < 0.05 & dat.response(:, 1) < 0));
        repSignfic = length(find(dat.pvalue < 0.05 & dat.response(:, 1) > 0));
end

% save for group plots
savefast(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator), 'dat');
close all;
end