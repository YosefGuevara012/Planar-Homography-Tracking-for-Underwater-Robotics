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


function [Hnew, WarpedImage, norm_x, norm_res, iter, last_non_zero_x, last_non_zero_res] = TrackImageSL3(ReferenceImage, CurrentImage, H, tracking_param)

global DEBUG_LEVEL_2;
if(DEBUG_LEVEL_2)
	disp('TrackingSL3');
	keyboard;
end;

% Initialise Homography
Hnew = H;
residues = 0;
iter = 1;
x = 10000000; % So that initially the norm(x) in the while loop is large


norm_x = zeros(1,tracking_param.max_iter);
norm_res = zeros(1,tracking_param.max_iter);

bool_norm = 1;
% Iterative minimization
%while(YOUR STOPPING CRITERION HERE)%
% (norm(x)> 0.01) & (iter <= 100)

% (norm(x)> tracking_param.max_x) & (iter <= tracking_param.max_iter) & (norm(rescale(residues))< tracking_param.max_err)
% (((norm(x)> tracking_param.max_x) || (norm(rescale(residues))< tracking_param.max_err)) && (iter <= tracking_param.max_iter))
% ((bool_norm == 1) && (norm(rescale(residues))< tracking_param.max_err) && (iter <= tracking_param.max_iter))
while(((norm(x)> tracking_param.max_x) || (norm(rescale(residues))< tracking_param.max_err)) && (iter <= tracking_param.max_iter)) 
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
%             set(fig, 'Position', [100, 100, 1280, 720]); % set figure size to 1280 x 720 pixels
			subplot(2,2,4); imagesc(WarpedImage.I); colormap(gray); title('Warped Image'); axis off;  
			W = zeros(ReferenceImage.sIv, ReferenceImage.sIu);
			W(WarpedImage.index) = weights;
			subplot(2,2,3); imagesc(W); colormap(gray); title('Weight Image'); axis off;
			subplot(2,2,2); image(abs(WarpedImage.I-double(ReferenceImage.I)).*WarpedImage.Mask.*W); 
			colormap(gray); 	title('Error Image'); 	axis off;  
			%keyboard;
            
            
		end; 
        
        norm_x(iter) = norm(x);
        norm_res(iter) = norm(rescale(residues));
        
        
		if(tracking_param.display)       
 
            fig3 = figure(3);
%             set(fig3, 'Position', [100, 100, 1280, 720]); % set figure size to 1280 x 720 pixels
            tiledlayout(2,1)
            
            
            %Add the new point to the animated line
            
            ax1 = nexttile;
            plot(ax1,norm_x(:),'r');
            grid on
            title('Euclidean norm of x')
            xlabel('Iteraction') 
            ylabel('norm(x)')
            
            ax2 = nexttile;
            plot(ax2,norm_res(:),'r');
            grid on
            title('Euclidean norm of the residues')
            xlabel('iteraction') 
            ylabel('norm(residues)')

            % Update the plot
            drawnow;
            %pause(0.5)
           
		end;
        
%         if (iter > 1)
%             if (norm_x(iter) > norm_x(iter-1))
%                 bool_norm = 0
%             end;
%             
%         end;
        
        
        iter = iter+1;

end;

if(tracking_param.display)

	% Transform polygon 
	WarpedImage.polygon = Hnew*ReferenceImage.polygon; 
	WarpedImage.polygon = WarpedImage.polygon./repmat(WarpedImage.polygon(3,:),3,1);

	norm(x);
	iter;

end

        last_non_zero_x = max(norm_x(find(norm_x, 1, 'last')));
        last_non_zero_res = max(norm_res(find(norm_res, 1, 'last')));


return

