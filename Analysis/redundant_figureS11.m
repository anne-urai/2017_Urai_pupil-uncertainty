function [] =  figureS11()
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

%% show two motionenergy filters to describe the filtering process
global mypath;

% get a random dataset
try
    load(sprintf('%s/Data/P17/Behav/P17_s1_2015-02-02_17-47-08.mat', mypath));
catch
    load('/Users/anne/Data/4Klaus/P7_s5_b10_2014-06-19_20-53-38.mat');
end
close all

% to avoid window UI bugs
display.dist        = window.dist;
display.res         = window.res;
display.width       = window.width;
display.frameRate   = window.frameRate;
display.center      = window.center;
clear sound audio

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute some general things for this file
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.ppd             = deg2pix(display, 1);  % pixel per degree
cfg.srange          = [-0.7 0.7];           % predetermine size the filter  space will have
cfg.frameRate       = display.frameRate;

% to avoid bugs with frameRates that are slightly below 60 Hz
if cfg.frameRate < 60 && round(cfg.frameRate)==60,
    cfg.frameRate = 60.1;
end

cfg.trange          = [0 0.2];              % temporal range
cfg.filtsize        = [ceil(diff(cfg.srange)*cfg.ppd) ...
    ceil(diff(cfg.srange)*cfg.ppd) ...
    ceil(diff(cfg.trange)*cfg.frameRate)];

% pad with some pixels for valid convolution
cfg.stimpad         = 50;
cfg.stimsize        = [2*dots.radius+1+cfg.stimpad 2*dots.radius+1+cfg.stimpad setup.nframes]; %  predetermine the size the cloud of dots will have
cfg.n_convolution   = cfg.stimsize + cfg.filtsize - 1;  % 'full'  size of the convolution
cfg.nextpow2        = 2.^nextpow2(cfg.n_convolution);
cfg.validsize       = cfg.stimsize - cfg.filtsize + 1;  % 'valid' size of convolution

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate filters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

direc                 = dots.direction;
counterdir            = dots.direction+180;
if counterdir > 360, counterdir = counterdir - 360; end
theta                 = [direc counterdir]; % only those two directions
theta = [90]; % up, for visualization

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. CREATE SPATIAL AND TEMPORAL FILTERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all parameters from Kiani et al., http://www.jneurosci.org/content/28/12/3017.full

% time and space axes
x = cfg.srange(1) : 1/cfg.ppd : cfg.srange(2);  % range in space, in degree, spaced by the size in degree of one pixel
y = cfg.srange(1) : 1/cfg.ppd : cfg.srange(2);  % range in space, in degree, spaced by the size in degree of one pixel
t = cfg.trange(1) : 1/cfg.frameRate: cfg.trange(2);   % range in time, in s, spaced by framerate

x = linspace(cfg.srange(1), cfg.srange(2), cfg.filtsize(1));
y = linspace(cfg.srange(1), cfg.srange(2), cfg.filtsize(2));
t = linspace(cfg.trange(1), cfg.trange(2), cfg.filtsize(3));

% important: use an uneven nr of points, easier for the later convolution
assert(mod(numel(x), 2) == 1, 'x should have an uneven nr of samples');
assert(mod(numel(y), 2) == 1, 'y should have an uneven nr of samples');
assert(mod(numel(t), 2) == 1, 't should have an uneven nr of samples');

% check if the filter has the size I thought it would get
assert(all([length(x) length(y) length(t)] == cfg.filtsize), 'filter does not have prespecified size');

% constants of the spatial filters
sc      = 0.35;
sg      = 0.05;

% create spatial mesh
[xmesh, ymesh] = meshgrid(x,y);

% constants of the temporal filters
k       = 60;
k = 85;
nslow   = 3;
nfast   = 5;

% TEMPORAL FUNCTIONS
% put them in the third dimension for the outer product
g1(1,1,:) = (k*t).^nslow .* exp(-k*t) .* (1/factorial(nslow) - ((k*t).^2) / factorial(nslow + 2));
g2(1,1,:) = (k*t).^nfast .* exp(-k*t) .* (1/factorial(nfast) - ((k*t).^2) / factorial(nfast + 2));

% rotate the meshgrid

yprime  = - xmesh .* sind(90) + ymesh .* cosd(90);
xprime  =   xmesh .* cosd(90) + ymesh .* sind(90);

% SPATIAL FUNCTIONS two fourth order cauchy functions
% transpose to match my unit circle directionality
alpha      = atand(xprime ./ sc);
f1rot      = (cosd(alpha).^4 .* cosd(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
f2rot      = (cosd(alpha).^4 .* sind(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';


%% make a nice looking schematic
figure;
colormap bone;
% start by showing the 2 temporal and 2 spatial filters

% temporal filters
subplot(441);
plot([t(1) t(end)], [0 0], 'k-'); hold on;
p = plot(t, squeeze(g1), t, squeeze(g2));
cols = linspecer(2);
p(1).Color = cols(2, :); p(2).Color = cols(1, :);

xlabel('Time (s)'); ylabel('Temporal filters');
ph2 = plot([0.2 0.21], 0.17*ones(2, 10), '.w');
lh = legend(ph2); % make t
lh.String = {'\color[rgb]{0.915294117647059,0.281568627450980,0.287843137254902} g_1', ...
    '\color[rgb]{0.363921568627451,0.575529411764706,0.748392156862745} g_2'};

box off; legend boxoff; axisNotSoTight; axis square;
set(gca, 'xtick', [0 0.1 0.2], 'ytick', [-0.1 0.2]); % offsetAxes;

% spatial filters
subplot(442);
plot([x(1) x(end)], [0 0], 'k-'); hold on;
plot([0 0], [-2.5 3], 'k-');
p = plot(x, sum(f1rot), x, sum(f2rot));
p(1).Color = cols(2, :); p(2).Color = cols(1, :);
xlabel('Space (degree)'); ylabel('Spatial filters');

ph2 = plot([0.7 0.7], 3*ones(2, 10), '.w');
lh = legend(ph2); % make t
lh.String = {'\color[rgb]{0.915294117647059,0.281568627450980,0.287843137254902} f_1', ...
    '\color[rgb]{0.363921568627451,0.575529411764706,0.748392156862745} f_2'};
legend boxoff;
 axis square; box off; axisNotSoTight;
set(gca, 'xtick', [-0.7 0 0.7], 'ytick', [-2 3]); %offsetAxes;

for thistheta = 90,
    
    % rotate the meshgrid
    yprime  = - xmesh .* sind(thistheta) + ymesh .* cosd(thistheta);
    xprime  =   xmesh .* cosd(thistheta) + ymesh .* sind(thistheta);
    
    % SPATIAL FUNCTIONS two fourth order cauchy functions
    % transpose to match my unit circle directionality
    alpha      = atand(xprime ./ sc);
    f1rot      = (cosd(alpha).^4 .* cosd(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
    f2rot      = (cosd(alpha).^4 .* sind(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
    
    % these two linear filters are in space-time quadrature
    basis1   = bsxfun(@times, f1rot, g1) ;
    basis2   = bsxfun(@times, f1rot, g2);
    basis3  = bsxfun(@times, f2rot, g1);
    basis4  = bsxfun(@times, f2rot, g2);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot info about the filters
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % their 4 combinations
    
    colormap bone;
    subplot(4,4,5);
    imagesc(x, t, squeeze(sum(basis1, 1))', [-.4 .4]);
    set(gca, 'ytick', [0 0.1 0.2]); axis square; box off;
    set(gca, 'yticklabel', [], 'xticklabel', []);
    
    subplot(446);
    imagesc(x, t, squeeze(sum(basis2, 1))', [-.4 .4]);
    set(gca, 'ytick', [0 0.1 0.2]); axis square; box off;
    set(gca, 'yticklabel', [], 'xticklabel', []);
    
    subplot(447);
    imagesc(x, t, squeeze(sum(basis3, 1))', [-.4 .4]);
    set(gca, 'ytick', [0 0.1 0.2]); axis square; box off;
    set(gca, 'yticklabel', [], 'xticklabel', []);
    
    subplot(448);
    imagesc(x, t, squeeze(sum(basis4, 1))', [-.4 .4]);
    set(gca, 'ytick', [0 0.1 0.2]); axis square; box off;
    set(gca, 'yticklabel', [], 'xticklabel', []);
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a separate filter for each direction
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 9;

for thistheta = theta,
    
    % rotate the meshgrid
    yprime  = - xmesh .* sind(thistheta) + ymesh .* cosd(thistheta);
    xprime  =   xmesh .* cosd(thistheta) + ymesh .* sind(thistheta);
    
    % SPATIAL FUNCTIONS two fourth order cauchy functions
    % transpose to match my unit circle directionality
    alpha      = atand(xprime ./ sc);
    f1rot      = (cosd(alpha).^4 .* cosd(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
    f2rot      = (cosd(alpha).^4 .* sind(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
    
    % these two linear filters are in space-time quadrature
    filt1   = bsxfun(@times, f1rot, g1) + bsxfun(@times, f2rot, g2);
    filt2   = bsxfun(@times, f2rot, g1) - bsxfun(@times, f1rot, g2);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot info about the filters
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % their 4 combinations
    subplot(4,4,cnt);
    imagesc(x, t, squeeze(sum(filt1, 1))', [-.4 .4]);
    xlabel('Space (degree)');
    if cnt == 5, ylabel('Time (s)'); end
    set(gca, 'ytick', [0 0.1 0.2]); axis square; box off;
    cnt  = cnt + 1;
    
    subplot(4,4,cnt);
    imagesc(x, t, squeeze(sum(filt2, 1))', [-.4 .4]);
    xlabel('Space (degree)');
    set(gca, 'ytick', [0 0.1 0.2]); axis square; box off;
    cnt = cnt + 1;
end


print(gcf, '-dpdf', sprintf('%s/Figures/FigureS11.pdf', mypath));

end% function end

function pix = deg2pix(display, ang)

% Converts visual angles in degrees to pixels.
%
% Inputs:
% display.dist (distance from screen (cm))
% display.width (width of screen (cm))
% display.res (number of pixels of display in horizontal direction)
% ang (visual angle)
% Warning: assumes isotropic (square) pixels
%
% Written 11/1/07 gmb zre
% Adapted by Anne Urai, August 2013

pixSize = display.width/display.res.width;   %cm/pix

sz = 2*display.dist*tan(pi*ang/(2*180));  %cm

pix = round(sz/pixSize);   %pix

return

end

