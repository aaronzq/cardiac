function [Pts_new, Pts_value] = coor_transform(Pts, M, l, lineV1, lineV2, lineV3, l_length)
%% Output the pts in new coordinate system defined by line l
% Output pts are in the same order as input pts

% Pts: [row, col, depth] in image coordinates;
% M: magnitude map
% l: center line in image coordinates; Note: it has to be the same
% downsampling(ds) rate as Pts, M
% lineV1, lineV2, lineV3: coor vectors of l
% voxelSize: voxelSize of Pts, M
       
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