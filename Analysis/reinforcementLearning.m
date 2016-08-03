function reinforcementLearning
% different reinforcement learning models

global mypath;
rng shuffle;
fits = table();

% function evaluation params
options                 = optimset('Display', 'off') ;
options.MaxFunEvals     = 5000000;
options.MaxIter         = 500000;
options.TolX            = 0.00000001;
options.TolFun          = 0.00000001;
options.Robust          = 'on';

% fit a number of different reinforcement learning models
for sj = 1:27,
    
    data  = readtable(sprintf('%s/Data/CSV/2ifc_data_sj%02d.csv', mypath, sj));
    
    % scale the pupil to be between 0 and 0.5 (just like uncertainty)
    data.decision_pupil = zscore(data.decision_pupil)/10 + 0.25;
    
    for v = 1:6,
        
        % 1. alpha 2. mu 3. sigma
        lowerbound   = [-1 -5 0];
        upperbound   = [1 5 10];
        firstguess   = [0 0 2.5];

        % minimise the cost function
        tic;
        [individualparams, NlogL, exitflag] = fminsearchbnd(@(individualparams) ...
            rlModel(data, individualparams, version), ...
            firstguess, lowerbound, upperbound, options);
        toc;
        
        if exitflag < 0, warning('exitflag indicates non-convergence'); end
        
        % save this to large table
        fits = cat(1, fits, ...
            cat(2, cell2table({sj, sprintf('model%d', v)}, 'variablenames', {'subjnr', 'model'}), ...
            array2table([individualparams NlogL], 'variablenames', {'alpha', 'mu1', 'sigma1', 'NlogL'})));
    end
end

% save
writetable(fits, sprintf('%s/Data/CSV/reinforcementLearning.csv', mypath));

end
