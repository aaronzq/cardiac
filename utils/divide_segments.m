function M_segment = divide_segments(M, line, ref_pts, voxel_size)
    
    [v1,v2,v3,~] = axis_on_line(line, ref_pts);    
    temp = sqrt(sum ( (line(2:end,:)-line(1:end-1,:)).^2 , 2));
    line_length = zeros(size(line,1),1);
    for i=2:size(line_length)
        line_length(i) = sum(temp(1:i-1));
    end
    
    [M_r, M_c, M_d] = ind2sub(size(M), find(M>0));
    [M_pts_new, ~] = coor_transform([M_r, M_c, M_d]*diag(voxel_size), M(M>0), line, v1, v2, v3, line_length);
    M_segment = zeros(size(M));
    for p = 1 : size(M_pts_new,1)
        if M_pts_new(p,2) < 0.5 * line_length(end)
            if M_pts_new(p,1) > 0
                M_segment(M_r(p),M_c(p),M_d(p)) = 1;
            else
                M_segment(M_r(p),M_c(p),M_d(p)) = 2;
            end            
        else
            if M_pts_new(p,1) > 0
                M_segment(M_r(p),M_c(p),M_d(p)) = 4;
            else
                M_segment(M_r(p),M_c(p),M_d(p)) = 3;
            end
        end   
    end
end