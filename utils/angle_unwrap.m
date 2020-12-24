function angle_new = angle_unwrap(angle)
%% angle: 0-360 degree 

    adjustment = zeros(size(angle));
    temp1 = 0;
    temp2 = 0;
    for i = 2:length(angle)
        if angle(i)-angle(i-1) < -300
            temp1 = temp1 + 360;
            if temp2 > 0
                adjustment = adjustment - 360;
            end
            temp2 = temp2 + 1;
        elseif angle(i)-angle(i-1) > 300
            temp1 = temp1 - 360;
            if temp2 < 0
                adjustment = adjustment + 360;
            end
            temp2 = temp2 - 1;
        else
            
        end
        adjustment(i) = adjustment(i) + temp1;
    end
    
    angle_new = angle + adjustment;

end