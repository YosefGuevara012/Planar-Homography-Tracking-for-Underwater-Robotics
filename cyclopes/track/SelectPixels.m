%===================================================================================
%
% Copyright (C) 2010. All rights reserved.
%
% This sofware was developed at:
% CNRS/I3S
% 2000 Route des Lucioles
% 06903 Sophia Antipolis
%
% NAME: SelectPixels 
% METHOD: Choose best pixels by analysing gradient 
% PRE: Reference image 
%
% POST: - Index, Jacobian for selected pixels
% AUTHORS: Andrew Comport
% DATE: 1/3/2010
%	CONTACT: comport@i3s.unice.fr
%
%====================================================================================


function SelectPixels(ReferenceImage)	

global DEBUG_LEVEL_1;
%if(DEBUG_LEVEL_1)
	disp('SelectPixels');
	keyboard;
%end;


	[U, S, V] = svd(ReferenceImage.J,'econ');
	[Sorted, Sorted_indexes] = sort(abs(U), 'descend');
	npixels = size(ReferenceImage.index,1);
	column_counter = ones(tracking_params.size_x,1);
	sub_index = [];
	new_index = [];

	chosen_table = zeros(size(ReferenceImage.J, 1));

	% Select pixels one at a time
	for(i=1:tracking_params.npixel_select)
		% Select max from each column one at a time
		for(j=1:tracking_params.size_x)
			% Check number of pixels i doesnt pass the required quantity
			if(i==tracking_params.npixel_select)
				j = tracking_params.npixel_select;
			else
			
				% Search for next available pixel in current column
				while(chosen_table(column_counter(j))~=0 & column_counter(j)~=npixels)
					column_counter(j)=column_counter(j)+1;
				end;

				% If a pixel is available
				if(column_counter(j)~=npixels)
					sub_index(i,1) = Sorted_indexes(column_counter(j),j);
					new_index(i,1) = ReferenceImage.index(sub_index(i,1));
					chosen_table(index(i,1)) = 1;
					i=i+1; 
				end;
			end;
		end;
	end;

	ReferenceImage.index = new_index;
	ReferenceImage.visibility_index = find(ReferenceImage.index);
	ReferenceImage.J = ReferenceImage.J(sub_index);
	ReferenceImage.JI.u = ReferenceImage.JI.u(sub_index);
	ReferenceImage.JI.v = ReferenceImage.JI.v(sub_index);
	ReferenceImage.Jinv = ReferenceImage.Jinv(:,sub_index);
	ReferenceImage.Ju = ReferenceImage.Ju(sub_index,:);
	ReferenceImage.Jv = ReferenceImage.Jv(sub_index,:);
	ReferenceImage.Jinv = pinv(ReferenceImage.J);


	%Mask = zeros(size(ReferenceImage.I));
	%Mask(new_index) = ones(size(new_index,2),1) ;
	%imagesc(Mask);


