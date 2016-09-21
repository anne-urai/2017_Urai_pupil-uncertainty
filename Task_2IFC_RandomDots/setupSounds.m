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

function sound = setupSounds(setup, audio)

%% auditory feedback
setup.feedbackcounterbalance = [1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 ...
   1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2  ];

if setup.participant > 0,
    switch setup.feedbackcounterbalance(setup.participant),
        case 1
            sound.feedback.correct      = 880; % 150 ms, 880 Hz
            sound.feedback.incorrect    = 200; % 150 ms, 200 Hz
            sound.stimonset             = 440; % 50 ms, 440 Hz
        case 2
            sound.feedback.correct      = 200; % 150  ms, 880 Hz
            sound.feedback.incorrect    = 880; % 150 ms, 200 Hz
            sound.stimonset             = 440; % 50 ms, 440 Hz
    end
else % in case running a test
    sound.feedback.correct      = 880; % 150  ms, 880 Hz
    sound.feedback.incorrect    = 200; % 150 ms, 200 Hz
    sound.stimonset             = 440; % 50 ms, 440 Hz
end

[sound.tonebuf, sound.tonepos] = CreateAudioBuffer(CreateTone(sound.stimonset, 0.050, audio.freq), ...
    CreateTone(sound.feedback.correct, 0.150, audio.freq), ...
    CreateTone(sound.feedback.incorrect ,0.150, audio.freq) , ...
    CreateTone(sound.feedback.incorrect, 0.5, audio.freq));

PsychPortAudio('FillBuffer', audio.h, sound.tonebuf);

end