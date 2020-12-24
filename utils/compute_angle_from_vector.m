function ang = compute_angle_from_vector(v_x, v_y, base)
    %% compute the angle (degree 0 - 360) from vector (2d)
    % v_x: vector component in x direction
    % v_y: vector component in y direction
    % base: base direction, rotate countclockwise --> increase angle


    if (v_x > 0 && v_y > 0)
        ang = atan(v_y/v_x) * 180 / pi;
    elseif(v_x > 0 && v_y < 0)
        ang = atan(v_y/v_x) * 180 / pi + 360;
    elseif(v_x < 0)
        ang = atan(v_y/v_x) * 180 / pi + 180;
    elseif(v_x == 0)
        ang = ((v_x >= 0) - 0.5)*2*90;       
    end
    
    ang = ang - base;
    if ang < 0
        ang = ang + 360;
    end
end