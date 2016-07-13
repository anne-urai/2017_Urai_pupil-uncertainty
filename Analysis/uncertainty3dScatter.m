%% 3d correlation
mypath = '/Users/anne/Data/pupilUncertainty_FigShare';

for sj = 1:27,
    data = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % add a measure of repeat
    repetition          = (diff(data.resp) == 0);
    data.repeat         = [repetition; NaN];
    data.evidence       = abs(data.motionstrength);
    
    % measure of how much the stimulus tells the subject to repeat from
    % their last choice to this one
    % important: the logistic slope for this measure of stimulus repetition
    % is steeper than for the one below!
    stimrepeat          = data.stim - circshift(data.resp, -1);
    stimrepeat          = double(stimrepeat == 0);
    stimrepeat(find(stimrepeat == 0)) = -1;
    stimrepeat          = stimrepeat .* circshift(abs(data.motionstrength), -1);
    data.stimrepeat     = stimrepeat;
    
    % remove the repeat measures for trials with non-consecutive trial nrs
    trlDif                      = [diff(data.trialnr); 0];
    removeTrls                  = false(size(trlDif));
    removeTrls(trlDif < 1)      = true;
    removeTrls(trlDif > 1)      = true;
    data.stimrepeat(removeTrls) = NaN;
    data.repeat(removeTrls)     = NaN;
    
    % =================================================== %
    % correlate evidence, pupil and repetition behaviour
    % separately for correct and error trials
    % =================================================== %
    
    cors = [0 1];
    for c = 1:2,
        thisdat   = data((data.correct == cors(c)), :);
        
        % 1. evidence and pupil
        mdl = fitglm(thisdat, 'evidence ~ 1 + decision_pupil');
        grandavg.evidencePupil(sj, c) = mdl.Coefficients.Estimate(2);
        
        % 2. pupil and repeat
        mdl = fitglm(thisdat, 'repeat ~ 1 + decision_pupil', ...
             'distr', 'binomial', 'link', 'logit');
        grandavg.pupilRepeat(sj, c) = mdl.Coefficients.Estimate(2);
        
        % 3. evidence and repeat
             mdl = fitglm(thisdat, 'repeat ~ 1 + evidence', ...
             'distr', 'binomial', 'link', 'logit');
        grandavg.evidenceRepeat(sj, c) = mdl.Coefficients.Estimate(2);
    end
end

%%
close all;
colors = cbrewer('qual', 'Set1', 3);
for c = 1:2,
    ax(c) = subplot(2,2,c);
    
    x = grandavg.evidencePupil(:, c);
    y = grandavg.pupilRepeat(:, c);
    z = grandavg.evidenceRepeat(:, c);
    
    % best fitting plane
    sf = fit([x, y],z,'poly11');
    p = plot(sf);
    grid on;
    set(p,'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3)
    hold on;
    scatter3(x,y,z, [], colors(c, :), 'filled');
    axis tight;
    view(3);
    
    xlabel('evidence - pupil'); 
    ylabel('pupil - repeat'); 
    zlabel('evidence - repeat');
    
    switch c
        case 1
            title('Error');
        case 2
            title('Correct');
    end
end
set(findall(gcf,'type','text'),'FontSize',10);
linkprop([ax(1) ax(2)], 'CameraPosition'); %h is the axes handle    

saveas(gcf,sprintf('%s/Figures/uncertainty3dScatter.fig', mypath),'fig');

%% plot in 2d
for c = 1:2,
    
    dat.evidencePupil   = grandavg.evidencePupil(:, c);
    dat.pupilRepeat     = grandavg.pupilRepeat(:, c);
    dat.evidenceRepeat  = grandavg.evidenceRepeat(:, c);
    figure;
    corrplot(dat, {'evidencePupil', 'pupilRepeat', 'evidenceRepeat'});
    
    switch c
        case 1
            suplabel('Errors', 't');
        case 2
            suplabel('Correct', 't');
    end
end

% now look 