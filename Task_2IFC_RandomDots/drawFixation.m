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

