%% 2-interval forced coice random dots
% Thresholding code at different coherence differences
%
% Anne Urai, February 2015
% To be used by OJay Medina
% -----------------------------------------------------------------

clear all; close all; clc;
addpath('D:\USERS\AnneUrai\Commitment\stats');
cd('D:\USERS\AnneUrai\OjayMedina\');
pwd

% general setup
setup.cancel        = false; % becomes true if escape is pressed, will abort experiment (but save data)
setup.Eye           = true; % true if using EyeLink
startatblock        = 1; % if anything went wrong during the session, can easily restart at a later block
setup.startatblock  = startatblock;

% ask for subject number, default: 0
setup.participant       = input('Participant number? ');
if isempty(setup.participant), setup.participant = 0; end

setup.session           = input('Session? ');
if isempty(setup.session), setup.session = 0; end

if startatblock > 1, %
    % load the file that was generated at block 1
    cd Data/
    allfiles = dir(fullfile('*.mat'));
    for file = 1:length(allfiles),
        % find the file that matches the participant and session nr
        if strncmp(sprintf('P%d_s%d_', setup.participant, setup.session), allfiles(file).name, 6);
            thisfile = allfiles(file).name;
        end
    end
    load(thisfile);
    disp('PREVIOUSLY CREATED FILE LOADED IN');
    setup.startatblock  = startatblock;
    setup.cancel        = false;
end
cd('D:\USERS\AnneUrai\OjayMedina\');

%% Setup the PsychToolbox
window.dist             = 50; % viewing distance in cm , 60 in EEG lab
window.width            = 40; % physical width of the screen in cm, 53.5 for BENQ in EEG lab
window.skipChecks       = 1; % set to 1 to skip VBL tests and avoid warnings
[window, audio]         = SetupPTB(window); %load all the psychtoolbox things
%ListenChar(2); % will stop keyboard output from going into the command window

Screen('TextSize', window.h, 15);
Screen('TextFont', window.h, 'Trebuchet');
Screen('TextColor', window.h, [255 255 255] );

if setup.startatblock == 1, %
    % create the design and load stimuli
    DrawFormattedText(window.h, 'Loading...',  'center', 'center');
    Screen('Flip', window.h);
    disp('configuring');
    [setup, dots, fix, results, sound, flip] = configuration_main(window, audio, setup);
else
    % fill the sound buffer again
    disp('keeping old configuration');
    PsychPortAudio('FillBuffer', audio.h, sound.tonebuf);
end

% display instructions
if setup.participant > 0,
    switch mod(setup.participant, 2);
        case 0
            switch setup.feedbackcounterbalance,
                case 1
                    DrawFormattedText(window.h, ['Klik met je linkerhand op Z voor zwakkere, \n \n' ...
                        'en met je rechterhand op M voor sterkere tweede beweging. \n \n' ...
                        'Je hoort een hoog piepje voor goede antwoorden, \n \n' ...
                        'en een laag piepje voor fouten. \n \n\n \n\n \n Succes!'],  'center', 'center');
                case 2
                    DrawFormattedText(window.h, ['Klik met je linkerhand op Z voor zwakkere, \n \n' ...
                        'en met je rechterhand op M voor sterkere tweede beweging. \n \n' ...
                        'Je hoort een laag piepje voor goede antwoorden, \n \n' ...
                        'en een hoog piepje voor fouten. \n \n\n \n\n \n Succes!'],  'center', 'center');
            end
        case 1
            switch setup.feedbackcounterbalance
                case 1
                    DrawFormattedText(window.h, ['Klik met je linkerhand op Z voor sterkere, \n \n' ...
                        'en met je rechterhand op M voor zwakkere tweede beweging. \n \n' ...
                        'Je hoort een hoog piepje voor goede antwoorden, \n \n' ...
                        'en een laag piepje voor fouten. \n \n\n \n\n \n Succes!'],  'center', 'center');
                case 2
                    DrawFormattedText(window.h, ['Klik met je linkerhand op Z voor sterkere, \n \n' ...
                        'en met je rechterhand op M voor zwakkere tweede beweging. \n \n' ...
                        'Je hoort een laag piepje voor goede antwoorden, \n \n' ...
                        'en een hoog piepje voor fouten. \n \n\n \n\n \n Succes!'],  'center', 'center');
            end
    end
end
Screen('Flip', window.h);
WaitSecs(.1);
KbWait();
WaitSecs(.1);

%% Start looping through the blocks trials
for block = setup.startatblock:setup.nblocks,
    
    DrawFormattedText(window.h, 'Loading...',  'center', 'center');
    Screen('Flip', window.h);
    
    % preload all the dot coordinates - overwrite this block!
    coord.ref           = nan(1, setup.ntrials, setup.nframes, 2, dots.nDots);
    coord.stim          = nan(1, setup.ntrials, setup.nframes, 2, dots.nDots);
    
    for trial = 1:setup.ntrials,
        % preload all the dot coordinates before starting the trial
        coord.ref(1, trial, :, :, :)        = dots_refstim(setup, window, dots, block, trial);
        coord.stim(1, trial, :, :, :)       = dots_limitedlifetime(setup, window, dots, block, trial);
    end
    
    % save all the dot coordinates
    save(sprintf('Data/Dots_P%d_s%d_b%d_%s.mat', setup.participant, setup.session, block, datestr(now, 'yyyy-mm-dd_HH-MM-SS')), '-mat', '-v7.3');
    
    if setup.Eye == true,
        edfFile = ELconfig(window, setup, block);
    end
    
    %% start the loop over trials
    for trial = 1:setup.ntrials,
        
        if setup.Eye,
            Eyelink ('Message', 'blinkbreak_start');
        end
        if trial == 1, % draw new dots, otherwise keep the ones from the last trial
            window      = dots_noise_draw(window, dots);
            window      = drawFixation(window, fix, dots); %fixation
            breaktime   = Screen('Flip', window.h); % flip once, so stationary dots
        end
        WaitSecs(setup.ISI);
        if setup.Eye,
            Eyelink ('Message', 'blinkbreak_end');
        end
        
        % stimulus sequence onset
        % FIXATION
        TimingCnt = GetSecs + window.frameDur - window.slack;
        
        for frameNum = 1:ceil(setup.fixtime(block, trial)*window.frameRate),
            
            window      = dots_noise_draw(window, dots);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.fix.VBL(block, trial, frameNum), ...
                flip.fix.StimOns(block, trial, frameNum), ...
                flip.fix.FlipTS(block, trial, frameNum), ...
                flip.fix.Missed(block, trial, frameNum), ...
                flip.fix.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.fix.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            
            if setup.Eye == true && frameNum == 1,
                Eyelink ('Message', sprintf('block%d_trial%d_fix', block, trial));
                fprintf('\n block%d_trial%d_fix \n', block, trial)
            end
        end
        
        % play reference stimulus onset tone
        PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
        PsychPortAudio('Start', audio.h); % like flip
        
        % REFERECE STIMULUS with 70% coherence
        for frameNum = 1:setup.nframes,
            
            window      = drawAllDots(window, dots, 1, trial, coord.ref, frameNum);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.refstim.VBL(block, trial, frameNum), ...
                flip.refstim.StimOns(block, trial, frameNum), ...
                flip.refstim.FlipTS(block, trial, frameNum), ...
                flip.refstim.Missed(block, trial, frameNum), ...
                flip.refstim.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.refstim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            
            if setup.Eye == true && frameNum == 1,
                Eyelink ('Message', sprintf('block%d_trial%d_ref70_dir%d', block, trial, dots.direction));
                fprintf('\n block%d_trial%d_ref \n', block, trial);
            end
        end
        
        % INTERVAL
        for frameNum = 1:ceil(setup.intervaltime(block, trial)*window.frameRate),
            
            window      = dots_noise_draw(window, dots);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.interval.VBL(block, trial, frameNum), ...
                flip.interval.StimOns(block, trial, frameNum), ...
                flip.interval.FlipTS(block, trial, frameNum), ...
                flip.interval.Missed(block, trial, frameNum), ...
                flip.interval.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.interval.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            
        end
        
        % play test stimulus onset tone
        PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
        PsychPortAudio('Start', audio.h); %like flip
        
        % TEST STIMULUS
        for frameNum = 1:setup.nframes,
            window      = drawAllDots(window, dots, 1, trial, coord.stim, frameNum);
            window      = drawFixation(window, fix, dots);
            
            [flip.stim.VBL(block, trial, frameNum), ...
                flip.stim.StimOns(block, trial, frameNum), ...
                flip.stim.FlipTS(block, trial, frameNum), ...
                flip.stim.Missed(block, trial, frameNum), ...
                flip.stim.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.stim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            
            if setup.Eye == true && frameNum == 1,
                Eyelink ('Message', sprintf('block%d_trial%d_stim_coh%.4f_dir%d_diff%d', ...
                    block, trial, dots.coherence(block, trial), dots.direction, setup.difficulty(block, trial)));
                fprintf('block%d_trial%d_stim_coh%.4f_dir%d_diff%d \n', ...
                    block, trial, dots.coherence(block, trial), dots.direction, find(round(abs(dots.coherence(block, trial)-0.7)*1000)== round(setup.cohlevels*1000)));
            end
        end
        
        %% RESPONSE
        frameNum = 1; keyIsDown = 0;
        while frameNum < setup.resptime*window.frameRate && ~keyIsDown,
            % when no response has been given, and the maximum response time hasnt been reached
            
            window      = dots_noise_draw(window, dots);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.resptime.VBL(block, trial, frameNum), ...
                flip.resptime.StimOns(block, trial, frameNum), ...
                flip.resptime.FlipTS(block, trial, frameNum), ...
                flip.resptime.Missed(block, trial, frameNum), ...
                flip.resptime.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.resptime.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            [keyIsDown, secs, keyCode]  = KbCheck();
            frameNum = frameNum + 1;
        end %button pressed
        
        results.resptime(block, trial)      = secs;
        
        if keyIsDown,
            
            results.RT(block, trial)            = results.resptime(block, trial) - flip.stim.VBL(block,trial, setup.nframes);
            results.press{block, trial}         = keyCode;
            
            try
                switch mod(setup.participant, 2),
                    case 0
                        switch KbName(keyCode),
                            case 'z', % left target 1, right target 2
                                results.response(block, trial) = -1;
                            case 'm',
                                results.response(block, trial) = 1;
                            case 'ESCAPE', % if escape is pressed, exit the experiment
                                setup.cancel = true;
                                results.response(block, trial) = NaN;
                            case 'esc', % if escape is pressed, exit the experiment
                                setup.cancel = true;
                                results.response(block, trial) = NaN;
                            otherwise % if any other key was pressed, fill in a NaN
                                results.response(block, trial) = NaN;
                        end
                    case 1  %Screen('DrawText',window.h,  'Press left for stronger and right for weaker motion',  window.center(1)*0.60, window.center(2)*0.55 , [255 255 255] );
                        switch KbName(keyCode),
                            case 'z', % left target 1, right target 2
                                results.response(block, trial) = 1;
                            case 'm',
                                results.response(block, trial) = -1;
                            case 'ESCAPE', % if escape is pressed, exit the experiment
                                setup.cancel = true;
                                results.response(block, trial) = NaN;
                            case 'esc', % if escape is pressed, exit the experiment
                                setup.cancel = true;
                                results.response(block, trial) = NaN;
                            otherwise % if any other key was pressed, fill in a NaN
                                results.response(block, trial) = NaN;
                        end
                end
            catch
                disp('did not recognize buttonpress');
            end
        end % buttonpress
        
        % code for correct responses
        if results.response(block, trial) == setup.increment(block,trial), %whether motion is stronger than 50% or not
            results.correct(block,trial) = true;
        elseif isnan(results.response(block, trial)),
            results.correct(block,trial) = NaN;
            results.RT(block,trial) = NaN; %set RT to NaN to easily filter out trials without a response
        else results.correct(block,trial) = false;
        end
        
        if setup.Eye == true,
            Eyelink ('Message', sprintf('block%d_trial%d_resp_key%d_correct%d', ...
                block, trial, results.response(block, trial), results.correct(block, trial)));
            fprintf('block%d_trial%d_resp_key%d_correct%d \n', ...
                block, trial, results.response(block, trial), results.correct(block, trial));
        end
        
        % pupil rebound
        % wait for the pupil to return to baseline, average 3s
        frameNum = 1;
        TimingCnt = GetSecs + window.frameDur - window.slack;
        
        while GetSecs < setup.pupilreboundtime(block, trial) + results.resptime(block, trial);
            window      = dots_noise_draw(window, dots);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.pupilrebound.VBL(block, trial, frameNum), ...
                flip.pupilrebound.StimOns(block, trial, frameNum), ...
                flip.pupilrebound.FlipTS(block, trial, frameNum), ...
                flip.pupilrebound.Missed(block, trial, frameNum), ...
                flip.pupilrebound.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.pupilrebound.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            frameNum = frameNum + 1;
        end
        
        %% FEEDBACK
        
        if results.correct(block,trial) == true, % correct
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
        elseif results.correct(block,trial) == false, % incorrect
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(3,1), sound.tonepos(3,2));
        elseif isnan(results.correct(block,trial)), % no response given, extra long tone
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(4,1), sound.tonepos(4,2));
        else %unrecognized response
            setup.cancel = true;
            warning('could not determine which feedback to give');
        end
        results.feedbackonset(block, trial) = GetSecs;
        results.soundstart.feedback(block, trial) = PsychPortAudio('Start', audio.h);
        
        if setup.Eye == true,
            Eyelink ('Message', sprintf('block%d_trial%d_feedback_correct%d_diff%d', ...
                block, trial, results.correct(block, trial),  find(round(abs(dots.coherence(block, trial)-0.7)*1000)== round(setup.cohlevels*1000))));
            fprintf('block%d_trial%d_feedback_correct%d_diff%d \n', ...
                block, trial, results.correct(block, trial),  find(round(abs(dots.coherence(block, trial)-0.7)*1000)== round(setup.cohlevels*1000)))
        end
        
        % capture the pupil response to feedback
        frameNum = 1;
        TimingCnt = GetSecs + window.frameDur - window.slack;
        
        while GetSecs < setup.pupilreboundtime2(block, trial) + results.feedbackonset(block, trial);
            window      = dots_noise_draw(window, dots);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.pupilrebound2.VBL(block, trial, frameNum), ...
                flip.pupilrebound2.StimOns(block, trial, frameNum), ...
                flip.pupilrebound2.FlipTS(block, trial, frameNum), ...
                flip.pupilrebound2.Missed(block, trial, frameNum), ...
                flip.pupilrebound2.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.pupilrebound2.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            frameNum = frameNum + 1;
        end
        % break out of all trials if ESC was pressed
        if setup.cancel,
            break
            warning('experiment was manually terminated');
        end
        
    end %end trial loop
    
    if block < setup.nblocks,
        DrawFormattedText(window.h, sprintf('Klaar met blok %d! \n \n Je had %.2f procent goed, \n \n en je antwoordde in gemiddeld %.2f seconden. \n \n\n \n Klik op enter om verder te gaan.', ...
            block, nanmean(results.correct(block,:))*100, nanmean(results.RT(block,:))), 'center', 'center');
    else % finish
        DrawFormattedText(window.h, sprintf('Klaar met blok %d! \n \n Je had %.2f procent goed, \n \n en je antwoordde in gemiddeld %.2f seconden. \n \n\n \n Einde van deze sessie!', ...
            block, nanmean(results.correct(block,:))*100, nanmean(results.RT(block,:))), 'center', 'center');
    end
    Screen('Flip', window.h);
    
    % also show this info in the command window
    fprintf('Finished block %d! \n \n You answered %.2f percent of trials correct, \n \n and your average reaction time was %.2f seconds. \n \n\n \n', ...
        block, nanmean(results.correct(block,:))*100, nanmean(results.RT(block,:)));
    
    % display accuracy over this block, per difficulty level
    fprintf('Correct per difficulty level: %.0f, %.0f, %.0f, %.0f, %.0f \n', ...
        100* nanmean(results.correct(block, abs(setup.coherence(block, :)-setup.cohlevels(1)) < 0.0001)), ...
        100* nanmean(results.correct(block, abs(setup.coherence(block, :)-setup.cohlevels(2)) < 0.0001)), ...
        100* nanmean(results.correct(block, abs(setup.coherence(block, :)-setup.cohlevels(3)) < 0.0001)), ...
        100* nanmean(results.correct(block, abs(setup.coherence(block, :)-setup.cohlevels(4)) < 0.0001)), ...
        100* nanmean(results.correct(block, abs(setup.coherence(block, :)-setup.cohlevels(5)) < 0.0001)));
    
    %% save the EL file for this block
    setup.datetime      = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    setup.eyefilename = sprintf('D:/USERS/AnneUrai/OjayMedina/Data/P%d_s%d_b%d_%s.edf', setup.participant, setup.session, block, setup.datetime);
    Eyelink('CloseFile');
    Eyelink('WaitForModeReady', 500);
    try
        status              = Eyelink('ReceiveFile',edfFile, setup.eyefilename); %this collects the file from the eyelink
        disp(status);
        disp(['File ' setup.eyefilename ' saved to disk']);
    catch
        warning(['File ' setup.eyefilename ' not saved to disk']);
    end

    % break out of all blocks if ESC was pressed
    if setup.cancel == true,
        break
        warning('experiment was manually terminated');
    end
    
end %end block loop

% wrap up and save
setup.datetime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% create subject specific file and save - add unique datestring to avoid any overwriting of files
setup.filename = sprintf('Data/P%d_s%d_%s.mat', setup.participant, setup.session, setup.datetime);
save(setup.filename, '-mat', 'setup', 'window', 'dots', 'fix', 'results', 'audio',  'sound', 'flip');
disp('SAVED FILE TO DISK'); disp(setup.filename);

% exit gracefully
disp('done!'); Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;
% ListenChar(1); % put keyboard output in the command window again

if block == setup.nblocks,
    % take all results and plot the full psychometric extravaganza
    [datapoints1, datapoints2, datapoints3, fit1, fit2, fit3, threshold] = FitThreshold_2IFC(setup.filename, 'fit', 'nobootstrap');
    
    % overwrite the perfrange file for the next time!
    save(sprintf('Data/P%d_threshold.mat', setup.participant), 'threshold');
    %save the fig for this session
    saveas(gcf, sprintf('Data/P%d_s%d_PsychFuncFit.pdf', setup.participant, setup.session), 'pdf');
end

%% show a measure of flip performance
figure;
for b = 1:setup.nblocks,
    for t = 1:setup.ntrials,
        try
            subplot(221);
            plot(diff(squeeze(flip.refstim.VBL(b,t,:))));
            if any(diff(squeeze(flip.refstim.VBL(b,t,:))) < 0.01 | diff(squeeze(flip.refstim.VBL(b,t,:))) > .02),
                fprintf('weird fliptime block %d, trial %d \n', b, t);
            end
            hold on;
        end
        try
            subplot(222);
            plot(diff(squeeze(flip.stim.VBL(b,t,:))));
            hold on;
        end
       % drawnow; %WaitSecs(0.05);
    end
end
subplot(221);
axis tight; title('int1'); xlabel('frames'); ylabel('frameDur');
axis tight; title('int2'); xlabel('frames'); ylabel('frameDur');