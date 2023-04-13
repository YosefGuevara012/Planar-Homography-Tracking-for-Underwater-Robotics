%===========================================================================
%
% Copyright (C) 2010. All rights reserved.
%
% This sofware was developed at:
% CNRS/I3S
% 2000 Route des Lucioles
% 06903 Sophia Antipolis
%
% NAME: Huber weighting function
% PRE: A robust measure of the scale of the error distribution (standard deviation) 
%			 The error vector
% POST: The weights corresponding to each residue 
% METHOD: See Robust Statistics, Peter. J. Huber, Wiley, 1981
% AUTHORS: Andrew Comport
% DATE: 1/1/2010
%	CONTACT: comport@i3s.unice.fr
%
%===========================================================================

function [weights, weights_index] = weightsHuber(scale, centered_residues)

% Huber threshold
c = 1.2107 * scale;

% Peicewise continuous function
center = find(centered_residues <  c);
tail = find(centered_residues >= c);

weights(center,1) = 1;
weights(tail,1) = c./centered_residues(tail);

weights_index = [center; tail];

return
