function [setup, dots, fix, results, sound, flip, trigger] = configuration_threshold_2IFC(window, audio, setup)

% general experimental design
setup.baselinecoh       = .7;
switch setup.perfrange
    case 'easy'
        setup.cohlevels         = [2.5, 5, 10, 20, 30] / 100;
    case 'medium'
        setup.cohlevels         = [1.25, 2.5, 5, 10, 30] / 100;
    case 'hard'
        setup.cohlevels         = [0.625, 1.25, 2.5, 5, 20] / 100;
end
setup.trialrep          = 100; %nr of trials per coherence level (per block)
setup.totalntrials      = round(setup.trialrep * length(setup.cohlevels));
setup.nblocks           = 10;
setup.ntrials           = setup.totalntrials / setup.nblocks;
assert(mod(setup.ntrials, 1)==0);

% timing
setup.fixtime           = (.5 + .5*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.viewingtime       = .75; % viewing duration in seconds (fixed in this script, or maximum viewing duration in RT paradigm
setup.nframes           = ceil(setup.viewingtime*window.frameRate); %number of frames the stimulus is displayed
setup.intervaltime      = (.3 + .4*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.resptime          = 3; % maximum time for response
setup.pupilreboundtime  = (1.5 + 1.5 * rand(setup.nblocks, setup.ntrials)); %wait average 4 seconds before displaying feedback
setup.pupilreboundtime2 = (2 + 1 * rand(setup.nblocks, setup.ntrials));
setup.ISI               = 1;

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
setup.increments         = [-1 1];
% use script from %Brooks, J.L. (2012). Counterbalancing for serial order carryover effects in 
% experimental condition orders. Psychological Methods.
% to counterbalance and randomize the order of responses

for b = 1:setup.nblocks,
    increm = carryoverCounterbalance(2,1,1+setup.ntrials/4,0); % 2 conditions, 1st order counterbalancing,
    increm(increm==2) = -1; % instead of condition nr 2, code for the -1 increment
    assert(length(increm)>=setup.ntrials, 'counterbalancing failed, not enough repetitions');
    
    setup.increment(b,:) = increm(1:setup.ntrials); % cut off the last repetitions if needed
end

setup.coherence          = Shuffle(repmat(setup.cohlevels', [ceil(setup.ntrials/length(setup.cohlevels)) setup.nblocks ]))';
setup.coherence          = setup.coherence(1:setup.nblocks, 1:setup.ntrials);

% for each trial, compute the actual coherence (from the baseline, individual threshold and the sign of the increment).
dots.coherence          = setup.baselinecoh + setup.increment.* setup.coherence;

%% sound
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
%flip.fix.StimOns            = flip.fix.VBL;
%flip.fix.FlipTS             = flip.fix.VBL;
%flip.fix.Missed             = flip.fix.VBL;
%flip.fix.beampos            = flip.fix.VBL;

flip.refstim.VBL            = nan(setup.nblocks, setup.ntrials, setup.nframes);
%flip.refstim.StimOns        = flip.refstim.VBL;
%flip.refstim.FlipTS         = flip.refstim.VBL;
%flip.refstim.Missed         = flip.refstim.VBL;
%flip.refstim.beampos        = flip.refstim.VBL;

flip.interval.VBL           = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.intervaltime))/window.frameDur));
%flip.interval.StimOns       = flip.interval.VBL;
%flip.interval.FlipTS        = flip.interval.VBL;
%flip.interval.Missed        = flip.interval.VBL;
%flip.interval.beampos       = flip.interval.VBL;

flip.stim.VBL               = nan(setup.nblocks, setup.ntrials, setup.nframes);
%flip.stim.StimOns           = flip.stim.VBL;
%flip.stim.FlipTS            = flip.stim.VBL;
%flip.stim.Missed            = flip.stim.VBL;
%flip.stim.beampos           = flip.stim.VBL;

flip.resptime.VBL           = nan(setup.nblocks, setup.ntrials, ceil(setup.resptime/window.frameDur));
%flip.resptime.StimOns       = flip.resptime.VBL;
%flip.resptime.FlipTS        = flip.resptime.VBL;
%flip.resptime.Missed        = flip.resptime.VBL;
%flip.resptime.beampos       = flip.resptime.VBL;

flip.pupilrebound1.VBL        = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.pupilreboundtime))/window.frameDur));
%flip.pupilrebound1.StimOns    = flip.pupilrebound1.VBL;
%flip.pupilrebound1.FlipTS     = flip.pupilrebound1.VBL;
%flip.pupilrebound1.Missed     = flip.pupilrebound1.VBL;
%flip.pupilrebound1.beampos    = flip.pupilrebound1.VBL;

flip.pupilrebound2.VBL        = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.pupilreboundtime2))/window.frameDur));
%flip.pupilrebound2.StimOns    = flip.pupilrebound2.VBL;
%flip.pupilrebound2.FlipTS     = flip.pupilrebound2.VBL;
%flip.pupilrebound2.Missed     = flip.pupilrebound2.VBL;
%flip.pupilrebound2.beampos    = flip.pupilrebound2.VBL;

end