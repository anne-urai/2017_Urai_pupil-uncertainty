function ang = pix2deg(window,pix)
%angle = pix2angle(window,pix)
%
%converts monitor pixels into degrees of visual angle.
%
%Inputs:
%window.dist (distance from screen (cm))
%window.width (width of screen (cm))
%window.resolution (number of pixels of window in horizontal direction)
%
%ang (visual angle)
%
%Warning: assumes isotropic (square) pixels

%Written 11/1/07 gmb zre

%Calculate pixel size
pixSize = window.width/window.res.width;   %cm/pix

sz = pix*pixSize;  %cm (duh)

ang = 2*180*atan(sz/(2*window(1).dist))/pi;


return

%test code

window.dist = 60; %cm
window.width = 44.5; %cm
window.resolution = [1680,1050];

pix = 100;

ang = pix2angle(window,pix)