function [window, audio] = SetupPTB(window)

% unify keycodes
KbName('UnifyKeyNames');

% skip PTB checks
if window.skipChecks,
    Screen('Preference', 'Verbosity', 0);
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'VisualDebugLevel', 0);
    % suppress warnings to the pput window
    Screen('Preference', 'SuppressAllWarnings', 1);
end

% % find the right screen
% screens = Screen('Screens');
% for s = 1:length(screens),
%    res = Screen('Resolution',screens(s));
%     if res.width == 1024 && res.height == 768 && res.hz == 60,
%         thisScreen = screens(s);
%     end
% end
% 
% res = Screen('Resolution',thisScreen);
% %assert((res.width == 1024 && res.height == 768 && res.hz == 60), 'DID NOT FIND CORRECT SCREEN');

thisScreen = 0;
% in lab D1.09, select the CRT monitor
window.screenNum = thisScreen;

% get the screen indices for the different colors
window.white    = WhiteIndex(window.screenNum);
window.black    = BlackIndex(window.screenNum);
window.gray     = round((window.white+window.black)/2); %rounding avoids problem with textures
window.bgColor  = window.black;
window.res      = Screen('Resolution', window.screenNum);

% Open the window
%[window.h, window.rect] =Screen('OpenWindow',window.screenNum,window.bgColor, [0 0 600 600]);
[window.h, window.rect] = Screen('OpenWindow',window.screenNum,window.bgColor);

% find out what happens without the blendfunction?
Screen(window.h,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Set the display parameters 'frameRate' and 'resolution'
window.frameDur = Screen('GetFlipInterval',window.h); %duration of one frame
window.frameRate = 1/window.frameDur; %Hz

[window.center(1), window.center(2)] = RectCenter(window.rect); % [window.rect(3)/2 window.rect(4)/2];
% include a slack for Flip 'when' management
window.slack    = window.frameDur / 3;

%% now the audio setup
InitializePsychSound(1);  % request low latency mode

audio           = [];
audio.freq      = 44100;
% open default soundport, in stereo (to match the sound matrix we create)
audio.h         = PsychPortAudio('Open', [], 1, [], audio.freq, 2, [], []);

HideCursor;
commandwindow;

disp('PTB setup complete');

% Show the subject some instructions
Screen('TextSize', window.h, 20);
Screen('DrawText',window.h, 'Loading the experiment.....', window.center(1)*0.60, window.center(2) , window.black );
Screen('Flip', window.h);

end