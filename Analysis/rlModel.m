function [NlogL, md] = rlModel(data, individualparams, version)
% ======================================= %
% applies reinforcement model fitting to  data in a number of different versions
%
%     v1: rpe (without belief state) updates perceptual noise
%     v2: rpe (without belief state) updates bound
%     v3: rpe (with belief state) updates perceptual noise
%     v4: rpe (with belief state) updates bound
%     v5: confidence updates perceptual noise
%     v6: confidence updates bound
%     v7: rpe (with pupil-based belief state) updates perceptual noise
%     v8: rpe (with pupil-based belief state) updates bound
%     v9: confidence updates bound
% 
% Anne Urai, 2016

% preallocate model struct
sz = nan(size(data, 1), 1);
md = struct('dv', sz, ...   % average DV, conditioned on the SJs choice
    'likelihood', sz, ...   % how often does the model choose as the sj?
    'p_strong', sz, ...     % chance the model would say 'strong' (0:1)
    'confidence', sz, ...   % expected value, 0.5:1
    'uncertainty', sz, ...  % uncertainty, 0:0.5
    'rpe', sz);             % rpe, reward - confidence

% only read this from the table once, saves time
md.stim        = data.motionstrength;   % including the external noise estimate
md.choice      = data.resp;             % md.choice(md.choice == -1) = 0;
md.reward      = data.correct;

if version > 6,
    md.pupil       = data.decision_pupil;
    % scale the pupil to be between 0 and 0.5 (just like uncertainty)
    md.pupil       = zscore(md.pupil)/10 + 0.25;
end

% learning rate will be fit
alpha      = individualparams(1);

% mean and sigma of internal noise will be updated iteratively
md.mu      = individualparams(2) * ones(size(md.stim));
md.sigma   = individualparams(3) * ones(size(md.stim));

% loop over trials
for t = 1:size(data, 1),
    
    % generate noise
    thisnoise = randn(1, 1000);
    
    % decision variable over a lot of different noisy samples
    dv             = md.stim(t) + (thisnoise * md.sigma(t));
    
    % let the model make a choice based on this information
    modelchoice    = sign(dv - md.mu(t));
    
    % find the trials in which the model matches the SJ's resp
    trls = find(modelchoice == md.choice(t));
    
    % how often does the model make the same choice as the subject?
    md.likelihood(t)     = length(trls) / length(modelchoice);
    
    % avoid breaking the solver (log(0) = Inf)
    if md.likelihood(t) == 0, md.likelihood(t) = 1 / size(thisnoise, 2); end
    
    % what's the chance the model would say 'strong'?
    modelchoice(modelchoice == -1) = 0;
    md.p_strong(t)      = mean(modelchoice);
    
    % do we use the belief state of the subject in their DV
    if version <= 2,
        md.dv(t)         = mean(dv);
    else
        md.dv(t)         = mean(dv(trls));
    end
    
    if version <= 6,
        % predict upcoming reward
        % md.confidence(t)     = 1 ./ (1 + exp(-slope .* (abs(md.dv(t) - bound)))); % Kahnt et al.
        md.confidence(t)     = 0.5 * (1 + erf(abs(md.dv(t) - md.mu(t)) ./ (md.sigma(t) * sqrt(2)))); % Lak, eq. 7
        
        % make sure the script doesn't crash on lapses
        if isnan(md.confidence(t));
            fprintf('.');
            md.confidence(t) = 0.5;
        end
        
    else
        % in v7 and v8, the pupil takes place of uncertainty
        md.confidence(t) = 1 - md.pupil(t);
    end
    
    % uncertainty = 1 - confidence
    md.uncertainty(t)   = 1 - md.confidence(t);
    
    % compute reward prediction error
    md.rpe(t) = md.reward(t) - md.confidence(t);
    
    % how and what does the model update?
    if t < size(data, 1),
        switch version
            case 1 % Kahnt style: rpe (without belief state) updates perceptual weight
                md.sigma(t+1) = md.sigma(t) + alpha * md.rpe(t);
                
            case 2 % classical RL: (without belief state) rpe updates value (=bound)
                md.mu(t+1) = md.mu(t) + alpha * md.rpe(t) * sign(md.choice(t));
                
            case 3 % rpe (with belief state) updates perceptual weight
                md.sigma(t+1) = md.sigma(t) + alpha * md.rpe(t);
                
            case 4 % Armin: (with belief state) rpe updates value (=bound)
                md.mu(t+1) = md.mu(t) + alpha * md.rpe(t) * sign(md.choice(t));
                
            case 5 % Kepecs + PL: confidence updates perceptual weight
                md.sigma(t+1) = md.sigma(t) + alpha * md.confidence(t);
                
            case 6 % confidence updates bound
                md.mu(t+1) = md.mu(t) + alpha * md.confidence(t) * sign(md.choice(t));
                
            case 7 % rpe pupil updates perceptual weight
                md.sigma(t+1) = md.sigma(t) + alpha * md.rpe(t);
                
            case 8 % rpe pupil  updates value (=bound)
                md.mu(t+1) = md.mu(t) + alpha * md.rpe(t) * sign(md.choice(t));
                
        end
    end % if t is not the last trial
    
end % trials

% this takes a long time
% modeldata = struct2table(md);

% compute NlogL on these data
NlogL       = -sum(log(md.likelihood));

% check that this will converge
assert(~isnan(NlogL));
assert(~isinf(abs(NlogL)));

% if any response wasnt modelled properly, we dont want these params
if any(md.confidence == 0.5),
    NlogL = NlogL + 100;
end

end % function
