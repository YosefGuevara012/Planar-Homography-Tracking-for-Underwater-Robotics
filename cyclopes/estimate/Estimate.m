%===========================================================================
%
% Copyright (C) 2010. All rights reserved.
%
% This sofware was developed at:
% CNRS/I3S
% 2000 Route des Lucioles
% 06903 Sophia Antipolis
%
% NAME: Estimate
% PRE: The Jacobian, the error vector and the tracking_param(m-estimator, method, noise threshold) 
% POST: The estimate and the weight vector 
% METHOD: Jacobian of I = Ic(w(P,x))
% AUTHORS: Andrew Comport
% DATE: 1/1/2010
%	CONTACT: comport@i3s.unice.fr
%
%===========================================================================

%function [x, weights] = Estimate(WarpedImage, ReferenceImage, residues, tracking_param)
function [x, weights] = Estimate(J, residues, tracking_param)

global DEBUG_LEVEL_3;
if(DEBUG_LEVEL_3)
	disp('Estimate');
	keyboard;
end;

% size of the data vector
n = size(residues,1);

% If M-estimator, need to recompute pseudo-inverse at each iteration
if(tracking_param.mestimator)

	% compute median centered residues and scale
	[scale, centered_residues] = mad(residues, tracking_param);

	% compte the weighting vector based on various function
 	% along with and index for errors completely rejected
	switch lower(tracking_param.robust_method)
	case 'huber'
		[weights, weights_index] = weightsHuber(scale, centered_residues);
	case 'tukey'
		[weights, weights_index] = weightsTukey(scale, centered_residues);
	end;

	W = repmat(weights,1,tracking_param.size_x);

	% Compute weighted pseudo-inverse and multiply by error vector
	%x = pinv(W(weights_index,:).*WarpedImage.J(weights_index,:))*(weights(weights_index).*residues(weights_index));
	x = pinv(W(weights_index,:).*J(weights_index,:))*(weights(weights_index).*residues(weights_index));
	
	% M-estimator with precomputed Jacobian...
	%x = ReferenceImage.Jinv(:,WarpedImage.visibility_index(weights_index))*(weights(weights_index).*residues(weights_index));

	%x = ReferenceImage.Jinv(:,weights_index).*residues(weights_index);
	%WJ = (W.*J);
	%x = pinv(WJ'*WJ)*WJ'*(w.*y);

else % Least squares estimation

	% Need to inverse Jacobian at each iteration if current or esm  
	if(tracking_param.estimation_method == 2 | tracking_param.estimation_method == 3)
		%x = pinv(WarpedImage.J)*residues;
		x = pinv(J)*residues;
	else
		 % Use precomputed pinv(J)
		 x = J*residues; % J is really Jinv in this case...
		 %x = Jinv(:,WarpedImage.visibility_index)*residues;
	end;

	weights = ones(size(n,1));

end;


return



