%===================================================================================
%
% Copyright (C) 2010. All rights reserved.
%
% This sofware was developed at:
% CNRS/I3S
% 2000 Route des Lucioles
% 06903 Sophia Antipolis
%
% NAME: InitTrackImageSL3 
% METHOD: Initialise reference image parameters for tracking
% PRE: Reference image 
%
% POST: - Mask for patch to be tracked
%				- Meshgrid for reference image
%				- Reference Image Jacobian and inverse
% AUTHORS: Andrew Comport
% DATE: 1/1/2010
%	CONTACT: comport@i3s.unice.fr
%
%====================================================================================


function ReferenceImage = InitTrackImageSL3(ReferenceImage, tracking_params);

global DEBUG_LEVEL_1;
if(DEBUG_LEVEL_1)
	disp('InitTrackSL3');
	keyboard;
end;

Jstack=[];

for i=1:tracking_params.ncolor

	% Compute constant reference Jacobian and keep the image gradient of the reference region 
	[ReferenceImage(i).J, ReferenceImage(i).JI] = JacobianImageSL3(ReferenceImage(i).I, ReferenceImage(i).P, ReferenceImage(i).index);

	% Compute geometric Jacobian (not computed separately in the previous function)
	[ReferenceImage(i).Ju ReferenceImage(i).Jv] = JacobianSL3(ReferenceImage(i).P, ReferenceImage(i).index);

	% Pixel selection
	if(tracking_params.pixel_selection == 1)
		SelectPixels(ReferenceImage(i));
	end;

	Jstack = [Jstack; ReferenceImage(i).J];

end;

% Pre-compute the pseduo inverse of the reference Jacobian
Jinv = pinv(Jstack);

start = 0;	
for i=1:tracking_params.ncolor

	finish = start + size(ReferenceImage(i).index,1);
	ReferenceImage(i).Jinv = Jinv(:,start+1:finish);
	start = finish;

end;

