function [ window ] = drawFixation( window, fix, dots )
% Fixation
% shape of the ABC (figure 1) in:
% Thaler, L., Sch?tz, a C., Goodale, M. a, & Gegenfurtner, K. R. (2013). 
% What is the best fixation target? The effect of target shape on stability
% of fixational eye movements. Vision research, 76, 31?42.

% cover the inner space for any stray dots
Screen('FillOval', window.h, window.black, [window.center(1)-dots.innerspace, ...
    window.center(2)- dots.innerspace*.8, window.center(1)+dots.innerspace*.8, window.center(2)+dots.innerspace*.8], ...
    dots.innerspace*.8); 

Screen('FillOval', window.h, fix.color, [window.center(1)-fix.circlesize/2, window.center(2)- fix.circlesize/2, window.center(1)+fix.circlesize/2, window.center(2)+fix.circlesize/2], fix.circlesize); 
Screen('DrawLine', window.h, window.black, window.center(1)-fix.circlesize/2, window.center(2), window.center(1)+fix.circlesize/2, window.center(2), fix.dotsize); 
Screen('DrawLine', window.h, window.black, window.center(1), window.center(2)-fix.circlesize/2, window.center(1), window.center(2)+fix.circlesize/2, fix.dotsize); 
% only the central dot changes color with feedback
Screen('FillOval', window.h, fix.color, [window.center(1)-fix.dotsize/2, window.center(2)- fix.dotsize/2, window.center(1)+fix.dotsize/2, window.center(2)+fix.dotsize/2], fix.dotsize); 

end

