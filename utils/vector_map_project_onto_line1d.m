function vels = vector_map_project_onto_line1d(U,V,W,line,voxel_size)
    
    line = unique(line,'row','stable');
    temp = line(2:end,:)-line(1:end-1,:);
    line_dir = [temp;0,0,0] + [0,0,0;temp];
    line_dir = line_dir ./ vecnorm(line_dir,2,2);
    
    [vecs_corr1, vecs_corr2, vecs_corr3] = ind2sub(size(U), find(U~=0 | V~=0 | W~=0));
    vecs_corr1 = vecs_corr1 * voxel_size(1);
    vecs_corr2 = vecs_corr2 * voxel_size(2);
    vecs_corr3 = vecs_corr3 * voxel_size(3);
    vec1 = U(U~=0 | V~=0 | W~=0);
    vec2 = V(U~=0 | V~=0 | W~=0);
    vec3 = W(U~=0 | V~=0 | W~=0);
    
    vels = zeros(size(line,1),1);
    for i=1:size(vecs_corr1,1)
        
        dist = vecnorm([line(:,1)-vecs_corr1(i),line(:,2)-vecs_corr2(i),line(:,3)-vecs_corr3(i)], 2, 2);
        [~, ind_sort] = sort(dist);
        min_ind = ind_sort(1);
        if vels(min_ind) == 0
            vels(min_ind) = line_dir(min_ind,:)*[vec2(i), vec1(i), vec3(i)]';  % row, col, depth
        else
            vels(min_ind) = 0.5*( vels(min_ind) + line_dir(min_ind,:)*[vec2(i), vec1(i), vec3(i)]');
        end
    end
    
    
end