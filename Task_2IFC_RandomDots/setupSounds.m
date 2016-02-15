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