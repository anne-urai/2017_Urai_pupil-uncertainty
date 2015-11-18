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

