function [trl, event] = trialfun_pupil(cfg)
% This code reproduces the analyses in the paper
% Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven
% by decision uncertainty and alters serial choice bias.
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai, 2016
% anne.urai@gmail.com

event   = cfg.event;
value   = {event(find(~cellfun(@isempty,strfind({event.value},'MSG')))).value};
sample  = [event(find(~cellfun(@isempty,strfind({event.value},'MSG')))).sample];

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * cfg.fsample);
posttrig =  round(cfg.trialdef.post * cfg.fsample);

trl = [];
session = cfg.session;
for j = 1:length(value), % loop through the trials and create the trial matrix on each trl
    
    % check that this is really a fixation trigger
    if (j < length(value) && ~isempty(strfind(value{j}, 'blinkbreak_end')) && ~isempty(strfind(value{j+2}, 'ref'))) || ...
            (j < length(value) && ~isempty(strfind(value{j+1}, 'fix')) && ~isempty(strfind(value{j+2}, 'ref'))) || ...
            (j == length(value) && ~isempty(strfind(value{j}, 'blinkbreak_end'))),
        
        trlbegin = sample(j) + pretrig;
        %trlend   = sample(j) + posttrig;
        offset   = pretrig;
        fixoffset = sample(j);
        
        % then find the trigger after that, corresponding to the reference
        if ~isempty(strfind(value{j+2}, 'ref')),
            refoffset = sample(j+2);
        else
            error('no refoffset sample found');
        end
        
        % find the trial nr and block nr, scan message
        scandat =  sscanf(value{j+2}, 'MSG %*f block%d_trial%d_ref70_*');
        blockcnt = scandat(1); trlcnt = scandat(2);
        
        % stimulus start
        if ~isempty(strfind(value{j+3}, 'stim_coh')),
            stimoffset = sample(j+3);
        else
            error('no stimoffset sample found');
        end
        
        % decode stimulus type
        coh =  sscanf(value{j+3}, ...
            'MSG %*f block%*d_trial%*d_stim_coh%f_dir%*d_diff%*d');
        if coh < .7,
            stimtype = -1; % weaker
        elseif coh > .7,
            stimtype = 1; % stronger
        elseif coh == .7,
            stimtype = 0; % no difference
        end
        
        % decode stim strength
        stimstrength = single(abs(.7-coh));
        if stimstrength >= 0.019 && stimstrength <= 0.03,
            stimstrength = 0.025;
        end
        
        % difficulty level
        diff =  sscanf(value{j+3}, ...
            'MSG %*f block%*d_trial%*d_stim_coh%*f_dir%*d_diff%d');
        if isempty(diff), diff = NaN; end
        
        % response
        if ~isempty(strfind(value{j+4}, 'resp')),
            respoffset = sample(j+4);
        else
            error('no respoffset sample found');
        end
        
        resp = sscanf(value{j+4}, ...
            'MSG %*f block%*d_trial%*d_resp_key%f_correct%f');
        if numel(resp) == 2,
            resptype = resp(1); respcorrect = resp(2);
        elseif numel(resp) == 0,
            % some erroneous asc files without key trigger
            respcorrect = sscanf(value{j+4}, 'MSG %*f block%*d_trial%*d_resp_key%*c_correct%d');
            switch respcorrect,
                case 1
                    resptype = stimtype;
                case 0
                    resptype = -stimtype;
            end
        end
        
        % feedback
        if length(value) > j+5 && ~isempty(strfind(value{j+5}, 'feedback')),
            feedbackoffset = sample(j+5);
            
            % check feedback type
            if ~isempty(strfind(value{j+5}, 'feedback_correct1')), % correct
                feedbacktype = 1;
            elseif ~isempty(strfind(value{j+5}, 'feedback_correct0')), % error
                feedbacktype = 0;
            elseif ~isempty(strfind(value{j+5}, 'feedback_correctNaN')), % no response given
                continue
                warning('no response trial removed');
            end
            
        else
            warning('no feedbackoffset sample found');
            global mypath;
            
            % load in the behavioural file that corresponds
            behavfile = dir(sprintf('%s/Data/Raw/P%02d/Behav/P%d_s%d_*.mat', ...
                mypath, cfg.sj, cfg.sj, cfg.session));
            
            % if there are several files, continue until the right one is
            % found....
            filefound = false; cnt = 1;
            while ~filefound,
                load(sprintf('%s/Data/Raw/P%02d/Behav/%s', mypath, cfg.sj, behavfile(cnt).name));
                try
                    if  ~all(isnan(results.correct(blockcnt, :))),
                        filefound = true;
                    end
                end
                cnt = cnt + 1;
            end
            assert(filefound == 1, 'could not find the right file');
            
            % get the feedbackoffset (= pupilrebound) between resp and tone
            feedbackoffset   = setup.pupilreboundtime(blockcnt, trlcnt);
            feedbackoffset   = round(respoffset + feedbackoffset*1000);
            feedbacktype     = results.correct(blockcnt, trlcnt);
        end
        
        % fieldtrip allows variable trial length
        trlend = feedbackoffset + posttrig;
        
        % append all to mimic the MEG's trialinfo
        newtrl   = [trlbegin trlend offset ...
            fixoffset refoffset stimtype stimstrength diff ...
            stimoffset resptype respcorrect respoffset ...
            feedbacktype feedbackoffset trlcnt blockcnt session];
        trl      = [trl; newtrl];
    end
    
end

% check that the coherence levels are properly coded
cohs = unique(trl(:,7));
if numel(cohs) ~= 5,
    error('difficulty levels not properly coded');
end

% if no difficulty was indicated in the triggers, insert based on the
% stimulus strength
for d = 1:length(cohs),
    trl(find(trl(:,7)==cohs(d)),8) = d;
end

end



