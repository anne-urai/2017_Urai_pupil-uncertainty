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

%% Motion strength example
% shows a random dot stimulus with a direction specific to the subject
% each trial decreases in motion strength

clear all; close all; clc;
addpath('D:\USERS\AnneUrai\Commitment\stats');
cd('D:\USERS\AnneUrai\OjayMedina\\');

% general setup
setup.cancel            = false;

% ask for subject number, default: 0
setup.participant       = input('Participant number? ');
if isempty(setup.participant), setup.participant = 0;
end

%% Setup the PsychToolbox
window.dist             = 50; % viewing distance in cm , 60 in EEG lab
window.width            = 40; % physical width of the screen in cm, 53.5 for BENQ in EEG lab
window.skipChecks       = 1; % set to 1 to skip VBL tests and avoid warnings
[window, audio]         = SetupPTB(window); %load all the psychtoolbox things

%% CONFIGURATION

[setup, dots, fix, results, sound, flip, coord] = configuration_motionstrength_example(window, audio, setup);

%% INSTRUCTIONS
Screen('TextSize', window.h, 15);
Screen('TextFont', window.h, 'Trebuchet');
Screen('TextColor', window.h, [255 255 255] );

DrawFormattedText(window.h, ['Welcome bij het experiment! \n \n '...
    'Je begint met een serie van instructie-schermen. \n\n' ...
    'Lees alle instructies goed, \n\n en klik op de spatiebalk om door te gaan naar het volgende scherm. \n\n' ...
    'Je kunt altijd vragen stellen als iets onduidelijk is.'],  'center', 'center');
Screen('Flip', window.h);
WaitSecs(.1); KbWait(); WaitSecs(.1);

DrawFormattedText(window.h, ['Het experiment bestaat uit herhalingen van hetzelfde taakje. \n\n' ...
    'In het begin staat steeds een wolk van witte puntjes op het scherm. \n\n \n\n' ...
    'In het midden zie je een rood fixatiepunt: \n\n het is belangrijk dat je hier je ogen op laat rusten, \n\n' ...
    'en geen oogbewegingen maakt naar de witte puntjes.'], ...
    'center', 'center');
Screen('Flip', window.h);
WaitSecs(.1); KbWait(); WaitSecs(.1);

window      = dots_noise_draw(window, dots);
window      = drawFixation(window, fix, dots); %fixation
Screen('Flip', window.h); % flip once, so stationary dots
WaitSecs(.1); KbWait(); WaitSecs(.1);

% do some randomization
switch mod(setup.participant, 4),
    case 0 %        dots.direction = 45;
        direction = 'rechter onderkant';
    case 1 %        dots.direction = 135;
        direction = 'linker onderkant';
    case 2 %        dots.direction = 225;
        direction = 'linker bovenkant';
    case 3 %        dots.direction = 315;
        direction = 'rechter bovenkant';
end

DrawFormattedText(window.h, ['Als het taakje start, hoor je een kort piepje \n\n' ...
    'en beginnen de witte puntjes te bewegen \n\n richting de ' direction ' van het scherm.'], 'center', 'center' )
Screen('Flip', window.h);
WaitSecs(.1); KbWait(); WaitSecs(.1);

%% Start looping through the blocks trials
for block = 1:setup.nblocks,
    
    %% start the loop over trials
    for trial = 1:setup.ntrials,
        
        if trial == 1, % draw new dots, otherwise keep the ones from the last trial
            window      = drawAllDots(window, dots, block, trial, coord.fix, 1);
            window      = drawFixation(window, fix, dots); %fixation
            Screen('Flip', window.h); % flip once, so stationary dots
        end
        
        WaitSecs(.1);
        
        %% stimulus sequence onset
        % FIXATION
        
        TimingCnt = GetSecs + window.frameDur - window.slack;
        for frameNum = 1:ceil(setup.fixtime(block, trial)*window.frameRate),
            
            window      = drawAllDots(window, dots, block, trial, coord.fix, frameNum);
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.fix.VBL(block, trial, frameNum), ...
                flip.fix.StimOns(block, trial, frameNum), ...
                flip.fix.FlipTS(block, trial, frameNum), ...
                flip.fix.Missed(block, trial, frameNum), ...
                flip.fix.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.fix.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        end
        
        % play reference stimulus onset tone
        PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
        PsychPortAudio('Start', audio.h); %like flip
        
        % TEST STIMULUS
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
        
        % break out of all trials if ESC was pressed
        if setup.cancel,
            break
            warning('experiment was manually terminated');
        end
        
        if trial == 1,
            WaitSecs(1);
            DrawFormattedText(window.h,  ['Zojuist kon je de bewegingsrichting van de puntjes duidelijk zien. \n\n' ...
                'Dat betekent dat de beweging STERK is, \n\n en dat alle puntjes dezelfde kant op gaan. '] , 'center', 'center');
            Screen('Flip', window.h);
            WaitSecs(.1); KbWait(); WaitSecs(.1);
            
            DrawFormattedText(window.h, ['Het kan ook zijn dat de puntjes meer flikkeren, \n\n'...
                'en dat de richting van beweging moeilijker te zien is. \n\n' ...
                'Dit betekent dat de beweging ZWAK is. '] , 'center', 'center');
            Screen('Flip', window.h);
            WaitSecs(.1); KbWait(); WaitSecs(.1);
            
            DrawFormattedText(window.h,  ['In dit experiment bewegen de puntjes altijd dezelfde kant op.\n\n ' ...
                'De STERKTE van de beweging zal verschillen: \n\n dit is belangrijk voor het taakje dat je gaat doen. \n\n' ...
                '\n\n Je ziet nu een aantal voorbeelden, \n\n beginnend met hele sterke en eindigend met hele zwakke beweging.' ], 'center', 'center')
            Screen('Flip', window.h);
            WaitSecs(.1); KbWait(); WaitSecs(.1);
        end
        
    end %end trial loop
    
    Screen('DrawText',window.h, 'Klaar!', window.center(1)*0.60, window.center(2) , [255 255 255] );
    Screen('Flip', window.h);
    WaitSecs(.1); KbWait(); WaitSecs(.1);
    
    % break out of all trials if ESC was pressed
    if setup.cancel,
        break
        warning('experiment was manually terminated');
    end
    
end %end block loop

% exit gracefully
disp('done!'); Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;
