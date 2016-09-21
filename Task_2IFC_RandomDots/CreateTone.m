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

function [finaltone] = CreateTone(tonefrequency,toneduration,samplingrate,attenuation,noise)
%% creates a tone with a pre-beep to trigger the UKE channel

if nargin < 5 || isempty(noise)
    noise = 0;
end
if nargin < 4 || isempty(attenuation)
    attenuation = 0;
end
if nargin < 3
    error('Not enough input arguments.');
end

% PRE-TONE 18KHZ BURST
dt = 0.005; % dampening length (s)
db = 25; % dampening attenuation (dB)
dn = dt*samplingrate;
sd = dn/sqrt(-2*log(10^(-db/20)));
pretonefrequency = 18000;
pretoneduration = .005;

pretone = sin(2*pi*pretonefrequency*(0:pretoneduration*samplingrate)/samplingrate);
if noise > 0
    pretone = pretone+noise*randn(size(pretone));
end

pretone = pretone*10^(-attenuation/20);

pretone(1:dn) = pretone(1:dn).*exp(-((1:dn)-dn).^2/(2*sd^2));
pretone(end-dn+1:end) = pretone(end-dn+1:end).*exp(-((1:dn)-1).^2/(2*sd^2));
pretone = min(max(pretone,-1),+1);
pretone = repmat(pretone,2,1);

% REAL TONE
dt = 0.010; % dampening length (s)
db = 40; % dampening attenuation (dB)
dn = dt*samplingrate;
sd = dn/sqrt(-2*log(10^(-db/20)));

tone = sin(2*pi*tonefrequency*(0:toneduration*samplingrate)/samplingrate);
if noise > 0
    tone = tone+noise*randn(size(tone));
end
tone = tone*10^(-attenuation/20);
tone(1:dn) = tone(1:dn).*exp(-((1:dn)-dn).^2/(2*sd^2));
tone(end-dn+1:end) = tone(end-dn+1:end).*exp(-((1:dn)-1).^2/(2*sd^2));
tone = min(max(tone,-1),+1);
tone = repmat(tone,2,1);

finaltone = [pretone tone];
%plot(finaltone')

end