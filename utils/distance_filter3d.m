function mask2 = distance_filter3d(mask,line,th)
           
    c_line = unique(round(line),'row','stable');
    mask2 = zeros(size(mask));
    for i = 1:length(c_line)
        mask2(c_line(i,1), c_line(i,2), c_line(i,3)) = 1;
    end
    se = strel('sphere',th);
    mask2 = convn(mask2, se.Neighborhood, 'same');
    mask2 = mask .* (mask2>0);

end