function [Pts_new, Pts_value] = coor_transform(Pts, M, l, lineV1, lineV2, lineV3, l_length)
%% Output new coordinates of input points, in reference to the line l 
% Preferably all inpupts in physical coordinates, rather than image coordinates

% Pts: coordinates in the order of [row, col, depth]
% M: magnitude map
% l: center line in the order of [row, col, depth]
% lineV1, lineV2, lineV3: coor vectors of l
% l_length: the length of line l
       
    Pts_new = zeros(size(Pts));
    Pts_value = zeros(size(Pts,1),1);
    for i=1:size(Pts,1)
        dist = vecnorm( l - Pts(i,:), 2, 2);
        [~, ind_sort] = sort(dist);
        min_ind = ind_sort(1);
        
        Pts_new(i,1) = (Pts(i,:) - l(min_ind,:)) * lineV1(min_ind,:)'; 
        Pts_new(i,2) = l_length(min_ind);
        Pts_new(i,3) = (Pts(i,:) - l(min_ind,:)) * lineV3(min_ind,:)';
               
        Pts_value(i) = M(i);
    end

end