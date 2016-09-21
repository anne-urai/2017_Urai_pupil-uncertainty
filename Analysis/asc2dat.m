function [data, event, blinksmp, saccsmp] = asc2dat(asc)
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

% create event structure for messages
evcell = cell(length(asc.msg),1);
event = struct('type', evcell, 'sample', evcell, 'value', evcell, 'offset', evcell, 'duration', evcell );

for i=1:length(asc.msg),
    
    strtok = tokenize(asc.msg{i});
    event(i).type = strtok{3};
    str2double(strtok{2});
    
    % match the message to its sample
    smpstamp = dsearchn(asc.dat(1,:)', str2double(strtok{2}));
    % find closest sample index of trigger in ascii dat
    
    if ~isempty(smpstamp)
        event(i).sample = smpstamp(1);
    else % if no exact sample was found
        warning('no sample found');
    end
    event(i).value = asc.msg{i};
end

% make data struct
% important: match the right data chans to their corresponding labels...
data                = [];
data.label          = {'EyeH'; 'EyeV'; 'EyePupil'};
data.trial          = {asc.dat(2:4, :)};  %% !!!!!!!!! %% only take gaze and pupil
data.fsample        = asc.fsample;
data.time           = {0:1/data.fsample:length(asc.dat(1,:))/data.fsample-1/data.fsample};
data.sampleinfo     = [1 length(asc.dat(1,:))];

if data.fsample ~= 1000,
    warning('pupil not sampled with 1000Hz');
end

% parse blinks
blinktimes = cellfun(@regexp, asc.eblink, ...
    repmat({'\d*'}, length(asc.eblink), 1), repmat({'match'}, length(asc.eblink), 1), ...
    'UniformOutput', false); % parse blinktimes from ascdat
blinktimes2 = nan(length(blinktimes), 2);
for s = 1:length(blinktimes), a = blinktimes{s};
    for j = 1:2, blinktimes2(s, j) = str2double(a{j}); end
end
timestamps = asc.dat(1,:); % get the time info
try
    blinksmp = arrayfun(@(x) find(timestamps == x, 1,'first'), blinktimes2, 'UniformOutput', true ); %find sample indices of blinktimes in timestamps
catch
    blinksmp = arrayfun(@(x) dsearchn(timestamps', x), blinktimes2, 'UniformOutput', true ); %find sample indices of blinktimes in timestamps
end

% parse saccades
sacctimes = cellfun(@regexp, asc.esacc, ...
    repmat({'\d*'}, length(asc.esacc), 1), repmat({'match'}, length(asc.esacc), 1), ...
    'UniformOutput', false); % parse blinktimes from ascdat
sacctimes2 = nan(length(sacctimes), 2);
for s = 1:length(sacctimes), a = sacctimes{s};
    for j = 1:2,
        if str2double(a{j}) ~= 0,
            sacctimes2(s, j) = str2double(a{j});
        else
            sacctimes2(s, j) = str2double(a{j+1});
        end
    end
end
try
    saccsmp = arrayfun(@(x) find(timestamps == x, 1,'first'), sacctimes2, 'UniformOutput', true ); %find sample indices
catch
    saccsmp = arrayfun(@(x) dsearchn(timestamps', x), sacctimes2, 'UniformOutput', true );
end

end