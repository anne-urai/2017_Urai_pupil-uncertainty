function [buffer,loop] = CreateAudioBuffer(varargin)

buffer = cat(2,varargin{:});

loop = zeros(length(varargin),2);
loop(1,:) = [0,size(varargin{1},2)-1];
for i = 2:length(varargin)
    loop(i,:) = loop(i-1,2)+[1,size(varargin{i},2)];
end

end