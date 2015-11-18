function pix = deg2pix(window, ang)

% Converts visual angles in degrees to pixels.
%
% Inputs:
% window.dist (distance from screen (cm))
% window.width (width of screen (cm))
% window.res (number of pixels of display in horizontal direction)
% ang (visual angle)
% Warning: assumes isotropic (square) pixels
%
% Written 11/1/07 gmb zre
% Adapted by Anne Urai, August 2013

pixSize = window.width/window.res.width;   %cm/pix

sz = 2*window.dist*tan(pi*ang/(2*180));  %cm

pix = round(sz/pixSize);   %pix 

return



