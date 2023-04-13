%===========================================================================
%
% Copyright (C) 2010. All rights reserved.
%
% This sofware was developed at:
% CNRS/I3S
% 2000 Route des Lucioles
% 06903 Sophia Antipolis
%
% NAME: Median Absolute Deviation 
% PRE: A vector of residues, and the number of parameters to estimate 
% POST: The robust scale of the error distribution and
%				The median centered residual error
% METHOD: MAD 
% AUTHORS: Andrew Comport
% DATE: 1/1/2010
%	CONTACT: comport@i3s.unice.fr
%
%===========================================================================

function [scale, centered_residues] = mad(residues, tracking_param)

n = size(residues,1);

% "median centered" residual error
centered_residues = abs(residues-median(residues));

% robust standard deviation (MAD)
%scale = 1.4826 * (1+5/(n-tracking_param.size_x)) * median(residues);
scale = 1.4826 * median(centered_residues);

%  Set a minimum threshold for scale
%  (when scale reaches the level of noise in the measurements)
 if(scale < tracking_param.scale_threshold)
 	%scale = 1.4826 * (1+5/(n-tracking_param.size_x)) * tracking_param.scale_threshold;
	 scale = 1.4826 * tracking_param.scale_threshold;
 end;

return

