%===============================================================
%
% Copyright (C) 2010. All rights reserved.
%
% This sofware was developed at:
% CNRS/I3S
% 2000 Route des Lucioles
% 06903 Sophia Antipolis
%
% NAME: WarpSL3
% METHOD: Homographic warping of image points with normalisation
% PRE: Points to be warped, Homography.
% POSE: Return warped points in indexed vector form
% AUTHOR: Andrew Comport
%	CONTACT: comport@i3s.unice.fr
%
%===============================================================

function WarpedImage = WarpSL3(ReferenceImage, H);

global DEBUG_LEVEL_3;
if(DEBUG_LEVEL_3)
	disp('WarpSL3');
	keyboard;
end;


% Small u,v for indexed vector form. Captial U,V for Matrix/Image form.
u = (H(1,1)*ReferenceImage.P.U(ReferenceImage.index)+...
		 H(1,2)*ReferenceImage.P.V(ReferenceImage.index)+H(1,3)); 
v = (H(2,1)*ReferenceImage.P.U(ReferenceImage.index)+...
		 H(2,2)*ReferenceImage.P.V(ReferenceImage.index)+H(2,3));
w = (H(3,1)*ReferenceImage.P.U(ReferenceImage.index)+...
		 H(3,2)*ReferenceImage.P.V(ReferenceImage.index)+H(3,3));

% Normalise
u = u./w;
v = v./w;

% Determine pixels inside of image and update indexes
WarpedImage.visibility_index = find(u<ReferenceImage.sIu+1 & u>=1 & v<ReferenceImage.sIv+1 & v>=1);
WarpedImage.index = ReferenceImage.index(WarpedImage.visibility_index);

% Stored in full image size matrix
WarpedImage.P.U = zeros(size(ReferenceImage.I));
WarpedImage.P.V = zeros(size(ReferenceImage.I));
WarpedImage.P.U(WarpedImage.index) = u(WarpedImage.visibility_index);
WarpedImage.P.V(WarpedImage.index) = v(WarpedImage.visibility_index);

% Propagate and update Mask and indexes
WarpedImage.Mask = zeros(ReferenceImage.sIv,ReferenceImage.sIu);
WarpedImage.Mask(WarpedImage.index) = 1;

return
