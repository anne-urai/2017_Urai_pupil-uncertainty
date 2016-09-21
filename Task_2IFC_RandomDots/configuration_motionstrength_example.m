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

function [setup, dots, fix, results, sound, flip, coord] = configuration_motionstrength_example(window, audio, setup)

% general experimental design
setup.cohlevels         = [95 90 80 70 60 50 40 30 20 10 0] / 100;
setup.trialrep          = 1; %nr of trials per coherence level (per block)
setup.totalntrials      = round(setup.trialrep * length(setup.cohlevels));
setup.nblocks           = 1;
setup.ntrials           = setup.totalntrials / setup.nblocks;
assert(mod(setup.ntrials, 1)==0);

% timing
setup.fixtime           = (.5 + .5*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.viewingtime       = 2; % viewing duration in seconds (fixed in this script, or maximum viewing duration in RT paradigm
setup.nframes           = ceil(setup.viewingtime*window.frameRate); %number of frames the stimulus is displayed
setup.resptime          = 10; % maximum time for response
setup.ISI               = 1;

%% design
[dots, fix]             = setupDots(window, setup);

% do some randomization
switch mod(setup.participant, 4),
    case 0
        dots.direction = 45;
    case 1
        dots.direction = 135;
    case 2
        dots.direction = 225;
    case 3
        dots.direction = 315;
end

dots.coherence  = setup.cohlevels;

%% auditory feedback
sound = setupSounds(setup, audio);

%% preallocate results and stimuli structure
% preallocation is a good habit to make sure that Matlab knows how big your
% output structures will be. You might run into memory problems if your
% structures grow on each loop - Matlab will have to find a new chunk of
% memory each time which costs significant time.

results.response            = NaN(setup.nblocks,setup.ntrials);
results.correct             = NaN(setup.nblocks,setup.ntrials);
results.RT                  = NaN(setup.nblocks,setup.ntrials);

% preallocate a full flip structure to store the output of every dynamic flip
flip.fix.VBL                = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.fixtime))/window.frameDur));
flip.fix.StimOns            = flip.fix.VBL;
flip.fix.FlipTS             = flip.fix.VBL;
flip.fix.Missed             = flip.fix.VBL;
flip.fix.beampos            = flip.fix.VBL;

flip.stim.VBL               = nan(setup.nblocks, setup.ntrials, setup.nframes);
flip.stim.StimOns           = flip.stim.VBL;
flip.stim.FlipTS            = flip.stim.VBL;
flip.stim.Missed            = flip.stim.VBL;
flip.stim.beampos           = flip.stim.VBL;

flip.resptime.VBL           = nan(setup.nblocks, setup.ntrials, ceil(setup.resptime/window.frameDur));
flip.resptime.StimOns       = flip.resptime.VBL;
flip.resptime.FlipTS        = flip.resptime.VBL;
flip.resptime.Missed        = flip.resptime.VBL;
flip.resptime.beampos       = flip.resptime.VBL;


%% preload all the dot coordinates
coord.fix           = nan(setup.nblocks, setup.ntrials, ceil(max(setup.fixtime(:))*window.frameRate), 2, dots.nDots);
coord.stim          = nan(setup.nblocks, setup.ntrials, setup.nframes, 2, dots.nDots);

for block = 1:setup.nblocks,
    for trial = 1:setup.ntrials,
        
        % preload all the dot coordinates before starting the trial
        coord.fix(block, trial, :, :, :)        = dots_noise(dots, ceil(max(setup.fixtime(:))*window.frameRate));
        coord.stim(block, trial, :, :, :)       = dots_limitedlifetime(setup, window, dots, block, trial);
    end
end

end