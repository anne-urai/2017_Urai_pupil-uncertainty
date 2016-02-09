function a6_retrieveDataFromPython2(whichmodulator)
%% this script requires that the intertrial toolbox in Python has been run for all participants
% read in the python-generated mat files and do some plots!

global mypath;
close all; clc;
subjects = 1:27;
nlags = 7;

for sj = subjects,
    
    % ============================================ %
    % ======= model WITH pupil term =========== %
    % ============================================ %
    
    load(sprintf('%s/Data/serialmodel/2ifc_%s_sj%02d.txtresults.mat', mypath, whichmodulator, sj));
    load(sprintf('%s/Data/serialmodel/2ifc_%s_sj%02d.txtdata.mat', mypath, whichmodulator, sj));
    
    switch whichmodulator
        case 'plain'
            model_h_mod = model_w_hist;
    end
    
    % get the fitted weights
    hf0 = model_h_mod.hf0+1;
    model_h_mod.w = model_h_mod.w(:);
    hlen = size(h, 2);
    
    % for bootstrapped values, also multiply with the bootstrapped slope
    bootstrap_corr(:,1:size(bootstrap,2)-2)  = bsxfun(@times, bootstrap(:, 1:end-2), bootstrap(:, end-1));
    
    weights = {'response', 'stimulus', ...
        'pupil', 'response_pupil', 'stimulus_pupil', ...
        'rt', 'response_rt', 'stimulus_rt'};
    cnt = 0;
    for w = 1:length(weights),
        thisw = model_h_mod.w(hf0+hlen*cnt:hf0+hlen*(cnt+1)-1);
        
        % project back into lag space
        thisw       = h * thisw;
        dat.(weights{w})(sj, :) = thisw;
        
        % get errorbars
        thisboot    = bootstrap_corr(:, nlags*cnt+1:nlags*(cnt+1));
        alpha       = 1 - 0.68; % should cover 1 std of the distribution
        thiswci     = prctile(thisboot, [100*alpha/2,100*(1-alpha/2)])';
        
        %thiswci(:, 1) = thisw - thiswci(:, 1);
        %thiswci(:, 2) = thiswci(:, 2) - thisw; % relative error
        dat.([weights{w} 'CI'])(sj, :, :) = thiswci;
        cnt = cnt  + 1;
    end
    
    % change into correct and error from stimulus and
    weights = {'', '_pupil', '_rt'};
    for w = 1:length(weights),        
        dat.(['correct' weights{w}])(sj, :) = ...
            dat.(['stimulus' weights{w}])(sj, :) + dat.(['response' weights{w}])(sj, :);
        dat.(['incorrect' weights{w}])(sj, :) = ...
            -dat.(['stimulus' weights{w}])(sj, :) + dat.(['response' weights{w}])(sj, :);
        
        dat.(['correct' weights{w} 'CI'])(sj, :, :) = ...
            dat.(['stimulus' weights{w} 'CI'])(sj, :, :) + dat.(['response' weights{w} 'CI'])(sj, :, :);
        dat.(['incorrect' weights{w} 'CI'])(sj, :, :) = ...
            -dat.(['stimulus' weights{w} 'CI'])(sj, :, :) + dat.(['response' weights{w} 'CI'])(sj, :, :);
    end
end

% save for group plots
savefast(sprintf('%s/Data/GrandAverage/historyweights_%s.mat', mypath, whichmodulator), 'dat');
close all;
end