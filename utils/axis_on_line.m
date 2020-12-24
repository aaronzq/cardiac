function [v1, v2, v3, v11] = axis_on_line(line, ref_pts)

v11 = repmat([ref_pts(1),ref_pts(2),ref_pts(3)], size(line,1), 1) - line;
temp = line(2:end,:) - line(1:end-1,:);
temp = temp ./ sqrt(sum(temp.^2,2));
v2 = [temp; 0 0 0;] + [0 0 0; temp;];
v3 = cross(v11, v2);
v1 = cross(v2, v3);

v11 = v11 ./ sqrt(sum(v11.^2,2));
v2 = v2 ./ sqrt(sum(v2.^2,2));
v3 = v3 ./ sqrt(sum(v3.^2,2));
v1 = v1 ./ sqrt(sum(v1.^2,2));

end