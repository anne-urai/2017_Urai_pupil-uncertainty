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