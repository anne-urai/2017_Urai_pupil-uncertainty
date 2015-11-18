function [ window ] = drawAllDots(window, dots, block, trial, stimuli, frameNum)

% draws the dots on the screen for this frame
Screen('DrawDots',window.h,squeeze(stimuli(block, trial, frameNum, :, :)), ...
    dots.size, dots.color, window.center, 2);
end