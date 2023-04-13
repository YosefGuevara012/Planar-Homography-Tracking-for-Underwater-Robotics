%===========================================================================
%
% Copyright (C) 2010. All rights reserved.
%
% This sofware was developed at:
% CNRS/I3S
% 2000 Route des Lucioles
% 06903 Sophia Antipolis
%
% NAME: Beaton-Tukey weighting function
% PRE: A robust measure of the scale of the error distribution (standard deviation) 
%			 The error vector
% POST: The weights corresponding to each residue 
% METHOD: See Robust Statistics, Peter. J. Huber, Wiley, 1981
% AUTHORS: Andrew Comport
% DATE: 1/1/2010
%	CONTACT: comport@i3s.unice.fr
%
%===========================================================================

function [weights, weights_index] = weightsTukey(scale, residues)

% Beaton-Tukey threshold
c = 4.6851 * scale;

weights_index = find(residues <  c);
outliers = find(residues >= c);

weights(weights_index, 1) = (1-(residues(weights_index)./c).^2).^2;
weights(outliers, 1) = 0.0; % completely cut-off outliers

return
