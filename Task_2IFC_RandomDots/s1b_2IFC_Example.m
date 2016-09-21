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

%% 2-interval forced coice random dots
% example with very easy motion strength before subjects start thresholding

clear all; close all; clc;
%addpath('D:\USERS\AnneUrai\Commitment\stats');
%cd('D:\USERS\AnneUrai\OjayMedina\');
pwd

% general setup
setup.cancel        = false; % becomes true if escape is pressed, will abort experiment (but save data)

% ask for subject number, default: 0
setup.participant       = input('Participant number? ');
if isempty(setup.participant), setup.participant = 0; end

%% Setup the PsychToolbox
window.dist             = 40; % viewing distance in cm , 60 in EEG lab
window.width            = 30; % physical width of the screen in cm, 53.5 for BENQ in EEG lab
window.skipChecks       = 1; % set to 1 to skip VBL tests and avoid warnings
[window, audio]         = SetupPTB(window); %load all the psychtoolbox things

%% CONFIGURATION

[setup, dots, fix, results, sound, flip, coord] = configuration_2IFC_example(window, audio, setup);

Screen('TextSize', window.h, 15);
Screen('TextFont', window.h, 'Trebuchet');
Screen('TextColor', window.h, [255 255 255] );

DrawFormattedText(window.h, ['Nu je weet hoe sterke en zwakke beweging eruit ziet,  \n \n' ...
    'ga je twee bewegingen met elkaar vergelijken.  \n\n' ...
    '\n Achter elkaar hoor je twee keer een piepje, \n \n en zie je twee keer beweging op het scherm.  \n \n\n\n' ...
    'Jouw taak is om te beslissen of de beweging de tweede keer  \n \n STERKER of ZWAKKER was dan de eerste keer. '],  'center', 'center');
Screen('Flip', window.h);
WaitSecs(.1); KbWait(); WaitSecs(.1);

switch mod(setup.participant, 2);
    case 0
        DrawFormattedText(window.h, ['Na het eind van de tweede beweging\n \n ' ...
            'geef je een ZWAKKER antwoord door de toets Z in te drukken met je linkerhand \n \n' ...
            'en geef je een STERKER antwoord door de toets M in te drukken met je rechterhand. \n \n\n\n' ...
            'Na je antwoord krijg je feedback in de vorm van een piepje,\n \n ' ...
            'en geschreven feedback op het scherm. \n \n'],  'center', 'center');
    case 1
        DrawFormattedText(window.h, ['Na het eind van de tweede beweging\n \n ' ...
            'geef je een STERKER antwoord door de toets Z in te drukken met je linkerhand \n \n' ...
            'en geef je een ZWAKKER antwoord door de toets M in te drukken met je rechterhand. \n \n\n\n' ...
            'Na je antwoord krijg je feedback in de vorm van een piepje,\n \n ' ...
            'en geschreven feedback op het scherm. \n \n'],  'center', 'center');
end
Screen('Flip', window.h);
WaitSecs(.1); KbWait(); WaitSecs(.1);

%% Start looping through the blocks trials
for block = 1:setup.nblocks,
    
    %% start the loop over trials
    for trial = 1:setup.ntrials,
        
        Screen('Flip', window.h); % flip once, so stationary dots
        window      = drawFixation(window, fix, dots); % fixation
        Screen('Flip', window.h); % flip once, so stationary dots
        WaitSecs(1);
        
        % ISI blinkbreak
        window      = dots_noise_draw(window, dots);
        window      = drawFixation(window, fix, dots); %fixation
        %        Screen('Flip', window.h); % flip once, so stationary dots
        %        WaitSecs(2);
        
        %% stimulus sequence onset
        % FIXATION
        TimingCnt = GetSecs + window.frameDur - window.slack;
        % window      = drawAllDots(window, dots, block, trial, coord.fix, 1);
        for frameNum = 1:ceil(setup.fixtime(block, trial)*window.frameRate),
            
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.fix.VBL(block, trial, frameNum), ...
                flip.fix.StimOns(block, trial, frameNum), ...
                flip.fix.FlipTS(block, trial, frameNum), ...
                flip.fix.Missed(block, trial, frameNum), ...
                flip.fix.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt, 1);
            TimingCnt = flip.fix.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        end
        
        % play reference stimulus onset tone
        PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
        results.soundstart.ref(block, trial) = PsychPortAudio('Start', audio.h); %like flip
        
        % REFERECE STIMULUS with 70% coherence
        TimingCnt = flip.fix.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        for frameNum = 1:setup.nframes,
            
            window      = drawAllDots(window, dots, block, trial, coord.ref, frameNum);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.refstim.VBL(block, trial, frameNum), ...
                flip.refstim.StimOns(block, trial, frameNum), ...
                flip.refstim.FlipTS(block, trial, frameNum), ...
                flip.refstim.Missed(block, trial, frameNum), ...
                flip.refstim.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.refstim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        end
        % INTERVAL
        TimingCnt = flip.refstim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        window      = drawAllDots(window, dots, block, trial, coord.interval, 1);
        for frameNum = 1:ceil(setup.intervaltime(block, trial)*window.frameRate),
            
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.interval.VBL(block, trial, frameNum), ...
                flip.interval.StimOns(block, trial, frameNum), ...
                flip.interval.FlipTS(block, trial, frameNum), ...
                flip.interval.Missed(block, trial, frameNum), ...
                flip.interval.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt, 1);
            TimingCnt = flip.interval.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        end
        
        % play test stimulus onset tone
        PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
        results.soundstart.stim(block, trial) = PsychPortAudio('Start', audio.h); %like flip
        
        % TEST STIMULUS
        TimingCnt = flip.interval.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        for frameNum = 1:setup.nframes,
            window      = drawAllDots(window, dots, block, trial, coord.stim, frameNum);
            window      = drawFixation(window, fix, dots);
            
            [flip.stim.VBL(block, trial, frameNum), ...
                flip.stim.StimOns(block, trial, frameNum), ...
                flip.stim.FlipTS(block, trial, frameNum), ...
                flip.stim.Missed(block, trial, frameNum), ...
                flip.stim.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.stim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        end
        
        %% RESPONSE
        TimingCnt = flip.stim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        keyIsDown = 0; frameNum = 1;
        window      = drawAllDots(window, dots, block, trial, coord.resp, 1);
        while frameNum < setup.resptime*window.frameRate && ~keyIsDown,
            % when no response has been given, and the maximum response time hasnt been reached
            
            window      = drawFixation(window, fix, dots); % fixation
            [flip.resptime.VBL(block, trial, frameNum), ...
                flip.resptime.StimOns(block, trial, frameNum), ...
                flip.resptime.FlipTS(block, trial, frameNum), ...
                flip.resptime.Missed(block, trial, frameNum), ...
                flip.resptime.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt, 1);
            TimingCnt = flip.resptime.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            
            frameNum = frameNum + 1;
            [keyIsDown, secs, keyCode]  = KbCheck();
            
        end %button pressed
        
        results.resptime(block, trial)      = GetSecs();
        
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
                disp('buttonpress not recognized');
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
        
        % pupil rebound
        % wait for the pupil to return to baseline, average 3s
        frameNum = 1;
        TimingCnt = GetSecs + window.frameDur - window.slack;
        % window      = dots_noise_draw(window, dots);
        
        while GetSecs < setup.pupilreboundtime(block, trial) + results.resptime(block, trial);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.pupilrebound2.VBL(block, trial, frameNum), ...
                flip.pupilrebound2.StimOns(block, trial, frameNum), ...
                flip.pupilrebound2.FlipTS(block, trial, frameNum), ...
                flip.pupilrebound2.Missed(block, trial, frameNum), ...
                flip.pupilrebound2.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt, 1);
            TimingCnt = flip.pupilrebound2.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            frameNum = frameNum + 1;
        end
        
        %% FEEDBACK
        
        if results.correct(block,trial) == true, % correct
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
            %  Screen('DrawText', window.h, 'GOED', window.center(1)*0.95, window.center(2) , [255 255 255] );
        elseif results.correct(block,trial) == false, % incorrect
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(3,1), sound.tonepos(3,2));
            %   Screen('DrawText', window.h, 'FOUT', window.center(1)*0.95, window.center(2) , [255 255 255] );
        elseif isnan(results.correct(block,trial)), % no response given
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(4,1), sound.tonepos(4,2));
            %  Screen('DrawText', window.h, 'GEEN ANTWOORD', window.center(1)*0.95, window.center(2) , [255 255 255] );
        else %unrecognized response
            setup.cancel = true;
            warning('could not determine which feedback to give');
        end
        results.feedbackonset(block, trial)       = GetSecs;
        results.soundstart.feedback(block, trial) = PsychPortAudio('Start', audio.h);
        %Screen('Flip', window.h);
        %WaitSecs(0.1);
        %
        % wait for the pupil to return to baseline, average 3s
        frameNum = 1;
        TimingCnt = GetSecs + window.frameDur - window.slack;
        % window      = dots_noise_draw(window, dots);
        
        while GetSecs < setup.pupilreboundtime2(block, trial) + results.feedbackonset(block, trial);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.pupilrebound2.VBL(block, trial, frameNum), ...
                flip.pupilrebound2.StimOns(block, trial, frameNum), ...
                flip.pupilrebound2.FlipTS(block, trial, frameNum), ...
                flip.pupilrebound2.Missed(block, trial, frameNum), ...
                flip.pupilrebound2.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt, 1);
            TimingCnt = flip.pupilrebound2.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            frameNum = frameNum + 1;
        end
        
        % break out of all trials if ESC was pressed
        if setup.cancel,
            break
            warning('experiment was manually terminated');
        end
        
        if 0,
            if trial == 1,
                WaitSecs(.5);
                DrawFormattedText(window.h, ['Onthoud de volgende dingen: \n \n \n \n \n \n' ...
                    '1. \n \n De beweging is beide keren even lang, \n \n' ...
                    'met een korte pauze ertussen. \n \n' ...
                    'In deze pauze is er geen informatie over bewegingssterkte te zien.\n \n '],  'center', 'center');
                Screen('Flip', window.h);
                WaitSecs(.1); KbWait(); WaitSecs(.1);
                
                DrawFormattedText(window.h, ['2. \n \n De taak is het makkelijkst als je je ogen op het rode punt in het midden houdt,' ...
                    '\n \n en individuele witte stipjes NIET probeert te volgen met je ogen. \n \n' ...
                    'Probeer de het totaalplaatje van de stipjeswolk te zien. '],  'center', 'center');
                Screen('Flip', window.h);
                WaitSecs(.1); KbWait(); WaitSecs(.1);
                
                DrawFormattedText(window.h, ['3.\n \n Geef pas antwoord als de tweede beweging is afgelopen. \n \n' ...
                    'Vanaf dat moment heb je drie seconden om te antwoorden.  \n \n \n\n' ...
                    'Als je te vroeg of te laat op de knop klikt, \n \n wordt je antwoord fout gerekend.  \n \n \n\n' ...
                    'Geef maar een 1x antwoord, en wacht dan op de feedback.'],  'center', 'center');
                Screen('Flip', window.h);
                WaitSecs(.1); KbWait(); WaitSecs(.1);
                
                DrawFormattedText(window.h, ['4. \n \n Luister goed naar de feedback piepjes: \n \n' ...
                    'je ziet alleen tijdens deze oefening de feedback ook op het scherm geschreven.'],  'center', 'center');
                Screen('Flip', window.h);
                WaitSecs(.1); KbWait(); WaitSecs(.1);
                
            elseif trial == 5,
                DrawFormattedText(window.h, ['Zoals je merkt, blijven de puntjes na feedback flikkeren. \n \n' ...
                    'Het is belangrijk dat je je ogen nog open houdt. \n \n' ...
                    'Als de puntjes stil staan, heb je een knipperpauze\n \n voordat je verdergaat met de volgende herhaling.'],  'center', 'center');
                Screen('Flip', window.h);
                WaitSecs(.1); KbWait(); WaitSecs(.1);
            end
        end
        
    end %end trial loop
    Screen('Flip', window.h);
    
    DrawFormattedText(window.h, sprintf('Klaar! \n \n Je had %.2f procent goed, \n \n en je antwoordde in gemiddeld %.2f seconden. \n \n\n \n Klik om verder te gaan.', ...
        nanmean(results.correct(block,:))*100, nanmean(results.RT(block,:))), 'center', 'center');
    
    Screen('Flip', window.h);
    WaitSecs(.1); KbWait(); WaitSecs(.1);
    
    fprintf('Klaar! \n \n Je had %.2f procent goed, \n \n en je antwoordde in gemiddeld %.2f seconden. \n \n\n \n Klik om verder te gaan. \n\n', ...
        nanmean(results.correct(block,:))*100, nanmean(results.RT(block,:)))
    % break out of all blocks if ESC was pressed
    if setup.cancel == true,
        break
        warning('experiment was manually terminated');
    end
    
    WaitSecs(.1); KbWait(); WaitSecs(.1);
    
end %end block loop

% exit gracefully
disp('done!'); Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;
