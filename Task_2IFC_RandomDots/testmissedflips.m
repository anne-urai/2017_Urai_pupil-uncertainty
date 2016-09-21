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

function testmissedflips(window, flip, setup)

cutoff = window.frameDur*1.5;
differenceflips = [];
differenceflipsstim = [];

%%% better:
for b = 1:setup.nblocks,
    for t = 1:setup.ntrials, %first trial of each block is discarded
        if exist('flip'),
            % put the whole series of flips for this one trial back together
            
            if setup.thresholding,
                flips = [squeeze(flip.fix.VBL(b,t, :))' squeeze(flip.stim.VBL(b,t, :))'];
            else
                flips = [squeeze(flip.fix.VBL(b,t, :))' squeeze(flip.refstim.VBL(b,t, :))' squeeze(flip.interval.VBL(b,t, :))' squeeze(flip.stim.VBL(b,t, :))'];
                
            end
            %flips = [ ];  squeeze(flip.pupilrebound2.VBL(b,t, :))'  squeeze(flip.resptime.VBL(b,t, :))'
            flips(isnan(flips)) = [];
            
            if ~isempty(find(diff(flips)>cutoff)); % the difference between flips on each trial
                disp(sprintf('FRAMES missed block %d trial %d', b, t));
            elseif ~isempty(find(diff(flips)<-cutoff)),
                disp(sprintf('FRAMES missed block %d trial %d', b, t));
                
            end
            
            differenceflips = [differenceflips diff(flips)];
            % differenceflipsstim = [differenceflipsstim diff(squeeze(flip.stim.VBL(b,t, :))')];
            
        end
    end
end % block

figure;
plot(differenceflips, '.b', 'MarkerSize', 20);
%plot(differenceflipsstim, '.r');

% title('Screen Flip performance');
xlabel('TRIALS'); ylabel('Time between Flips (MS)'); %ylim([0 0.02]);
set(gca, 'XTickLabel', get(gca, 'XTick')*window.frameDur);
title('PSYCHTOOLBOX FLIP TIMING TEST');

end

