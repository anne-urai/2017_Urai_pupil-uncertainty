function [setup, dots, fix, results, sound, flip, coord] = configuration_2IFC_example(window, audio, setup)

% general experimental design
setup.baselinecoh       = .7; %70% baseline coherence, great for MEG
setup.trialrep          = 20; %600; %nr of trials over the whole session - 600 for MEG, more for training?
setup.nblocks           = 1; % how many blocks? this way, one block takes about 7 min

setup.ntrials           = round(setup.trialrep/setup.nblocks); % total nr of trials per block
assert(mod(setup.ntrials, 1)==0);

% set a high motion difference for the example
setup.threshold         = .3;

% timing
setup.fixtime           = (.5 + .5*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.viewingtime       = 1; % viewing duration in seconds (fixed in this script, or maximum viewing duration in RT paradigm
setup.nframes           = ceil(setup.viewingtime*window.frameRate); %number of frames the stimulus is displayed
setup.intervaltime      = (.3 + .4*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.resptime          = 1; % maximum time for response
setup.pupilreboundtime  = (.5 + .5 * rand(setup.nblocks, setup.ntrials)); %1.5-2.5s after resp
setup.pupilreboundtime2 = (.5 + .5 * rand(setup.nblocks, setup.ntrials));% 2-2.5s after fb

%% design
[dots, fix]             = setupDots(window, setup);

% direction depends on the participant
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

% decrease or increase from baseline coherence
setup.increments         = [1 -1];

% parameters that vary on each trial
for b = 1:setup.nblocks,
    increm               = repmat(setup.increments, 1, setup.ntrials/length(setup.increments));
   % setup.increment(b,:) = increm(randperm(length(increm)));
    setup.increment(b,:) = increm;
end

% for each trial, compute the actual coherence (from the baseline, individual threshold and the sign of the increment).
dots.coherence      = setup.baselinecoh + setup.increment.* setup.threshold;

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

results.soundstart.ref      = NaN(setup.nblocks,setup.ntrials);
results.soundstart.stim     = NaN(setup.nblocks,setup.ntrials);
results.soundstart.feedback = NaN(setup.nblocks,setup.ntrials);

% preallocate a full flip structure to store the output of every dynamic flip
flip.fix.VBL                = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.fixtime))/window.frameDur));
flip.fix.StimOns            = flip.fix.VBL;
flip.fix.FlipTS             = flip.fix.VBL;
flip.fix.Missed             = flip.fix.VBL;
flip.fix.beampos            = flip.fix.VBL;

flip.refstim.VBL            = nan(setup.nblocks, setup.ntrials, setup.nframes);
flip.refstim.StimOns        = flip.refstim.VBL;
flip.refstim.FlipTS         = flip.refstim.VBL;
flip.refstim.Missed         = flip.refstim.VBL;
flip.refstim.beampos        = flip.refstim.VBL;

flip.interval.VBL           = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.intervaltime))/window.frameDur));
flip.interval.StimOns       = flip.interval.VBL;
flip.interval.FlipTS        = flip.interval.VBL;
flip.interval.Missed        = flip.interval.VBL;
flip.interval.beampos       = flip.interval.VBL;

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

flip.pupilrebound1.VBL        = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.pupilreboundtime))/window.frameDur));
flip.pupilrebound1.StimOns    = flip.pupilrebound1.VBL;
flip.pupilrebound1.FlipTS     = flip.pupilrebound1.VBL;
flip.pupilrebound1.Missed     = flip.pupilrebound1.VBL;
flip.pupilrebound1.beampos    = flip.pupilrebound1.VBL;

flip.pupilrebound2.VBL        = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.pupilreboundtime2))/window.frameDur));
flip.pupilrebound2.StimOns    = flip.pupilrebound2.VBL;
flip.pupilrebound2.FlipTS     = flip.pupilrebound2.VBL;
flip.pupilrebound2.Missed     = flip.pupilrebound2.VBL;
flip.pupilrebound2.beampos    = flip.pupilrebound2.VBL;

%% preload all the dot coordinates
coord.fix           = nan(setup.nblocks, setup.ntrials, ceil(max(setup.fixtime(:))*window.frameRate), 2, dots.nDots);
coord.ref           = nan(setup.nblocks, setup.ntrials, setup.nframes, 2, dots.nDots);
coord.interval      = nan(setup.nblocks, setup.ntrials, ceil(max(setup.intervaltime(:))*window.frameRate), 2, dots.nDots);
coord.stim          = nan(setup.nblocks, setup.ntrials, setup.nframes, 2, dots.nDots);
coord.resp          = nan(setup.nblocks, setup.ntrials, ceil(max(setup.resptime(:))*window.frameRate), 2, dots.nDots);
%coord.pupil1        = nan(setup.nblocks, setup.ntrials, ceil(max(setup.pupilreboundtime(:))*window.frameRate), 2, dots.nDots);
%coord.pupil2        = nan(setup.nblocks, setup.ntrials,ceil(max(setup.pupilreboundtime2(:))*window.frameRate), 2, dots.nDots);

for block = 1:setup.nblocks,
    for trial = 1:setup.ntrials,
        
        % preload all the dot coordinates before starting the trial
        coord.fix(block, trial, :, :, :)        = dots_noise(dots, ceil(max(setup.fixtime(:))*window.frameRate));
        coord.ref(block, trial, :, :, :)        = dots_refstim(setup, window, dots, block, trial);
        coord.interval(block, trial, :, :, :)   = dots_noise(dots, ceil(max(setup.intervaltime(:))*window.frameRate));
        coord.stim(block, trial, :, :, :)       = dots_limitedlifetime(setup, window, dots, block, trial);
        coord.resp(block, trial, :, :, :)       = dots_noise(dots, ceil(max(setup.resptime(:))*window.frameRate));
       % coord.pupil1(block, trial, :, :, :)     = dots_noise(dots, ceil(max(setup.pupilreboundtime(:))*window.frameRate));
       % coord.pupil2(block, trial, :, :, :)     = dots_noise(dots, ceil(max(setup.pupilreboundtime2(:))*window.frameRate));
    
    end
end

end