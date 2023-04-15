%===================================================================
%
% Copyright (C) 2010. All rights reserved.
%
% This sofware was developed at:
% CNRS/I3S
% 2000 Route des Lucioles
% 06903 Sophia Antipolis
%
% NAME: TrackImageSL3 - Non-linear iterative estimation of 3x3 planar 
%												homography between reference and current image
% PRE: - Reference image with I, P, JI, JP
%			 - Current image with I
%			 - Current estimate/prediction of Homographie H
% POST:  Hnew - New esimate of Homography
% AUTHOR: Andrew Comport
% DATE: 1/1/08
%	CONTACT: comport@i3s.unice.fr
%
%===================================================================


function [Hnew, WarpedImage] = TrackImageSL3(ReferenceImage, CurrentImage, H, tracking_param)

global DEBUG_LEVEL_2;
if(DEBUG_LEVEL_2)
	disp('TrackingSL3');
	keyboard;
end;

% Initialise Homography
Hnew = H;
residue = 0;
iter = 0;
x = 10000000; % So that initially the norm(x) in the while loop is large 


% Iterative minimization
%while(YOUR STOPPING CRITERION HERE)% 
while(1 & (norm(x)> 0.01) & (iter <= 200))
		% Current patch
    WarpedImage = WarpImageSL3(CurrentImage, ReferenceImage, Hnew);    

	  % Patch error/residue in vector form
    residues = double(ReferenceImage.I(WarpedImage.index)) - WarpedImage.I(WarpedImage.index);

		switch tracking_param.estimation_method

			case 1
				% If M-estimator is used then dont use pre computed Jacobian
				if(tracking_param.mestimator == 1)
					WarpedImage.J = ReferenceImage.J(WarpedImage.visibility_index,:); 
				end;
				%Pseudo inverse precalculated if no mestimator
				[x, weights] = Estimate(WarpedImage, ReferenceImage, residues, tracking_param);

		  case 2
				% Compute Current Jacobian
				[WarpedImage.J, WarpedImage.JI] = JacobianImageSL3(WarpedImage.I, WarpedImage.P, WarpedImage.index);
				[x, weights] = Estimate(WarpedImage, ReferenceImage, residues, tracking_param);

			case 3 % If second order minimisation (ESM) recompute current image gradient at each iteration
				% Note esm jacobian stored in WarpedImage
			  WarpedImage.J = JacobianImageESMSL3(WarpedImage, ReferenceImage); 
				[x, weights] = Estimate(WarpedImage, ReferenceImage, residues, tracking_param);

			otherwise
				error('Tracking estimation method does not exist');

		end;

		% Compute unknown parameters x 
    A = [x(5),x(3),x(1); x(4),-x(5)-x(6),x(2); x(7),x(8),x(6)]; 
    Hnew = Hnew*expm(A); % Compute the update of the Homography matrix using the exponential map
        
        
    
		if(tracking_param.display)
			fig = figure(1);
            set(fig, 'Position', [100, 100, 1280, 720]); % set figure size to 1280 x 720 pixels
			subplot(2,2,4); imagesc(WarpedImage.I); colormap(gray); title('Warped Image'); axis off;  
			W = zeros(ReferenceImage.sIv, ReferenceImage.sIu);
			W(WarpedImage.index) = weights;
			subplot(2,2,3); imagesc(W); colormap(gray); title('Weight Image'); axis off;
			subplot(2,2,2); image(abs(WarpedImage.I-double(ReferenceImage.I)).*WarpedImage.Mask.*W); 
			colormap(gray); 	title('Error Image'); 	axis off;  
			%keyboard;
            
            
		end; 

		iter = iter+1;
        
        
		if(tracking_param.display)       

			norm_x(iter) = norm(x);
            
            fig2 = figure(2);
            set(fig2, 'Position', [100, 100, 1280, 720]); % set figure size to 1280 x 720 pixels
            % Add the new point to the animated line
            plot(norm_x(:),'r');
            title('Euclidean norm of x')
            xlabel('iteraction') 
            ylabel('norm(x)') 

            % Update the plot
            drawnow;
            pause(0.3)
           
		end;
end;

if(tracking_param.display)

	% Transform polygon 
	WarpedImage.polygon = Hnew*ReferenceImage.polygon; 
	WarpedImage.polygon = WarpedImage.polygon./repmat(WarpedImage.polygon(3,:),3,1);

	norm(x)
	iter

end

return

