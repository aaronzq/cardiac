function plot_vec(pt, vec, length, c, head)

    x = [pt(1) pt(1)+length*vec(1)];
    y = [pt(2) pt(2)+length*vec(2)];
    z = [pt(3) pt(3)+length*vec(3)];
    line(x,y,z,'Color',c,'lineWidth',2);
    
    if head == 1 %xy plane   
        % plot an arrow head on xy plane    
        v1 = [vec(1),vec(2)]; v1 = v1 ./ norm(v1);
        v2 = [vec(2),-vec(1)]; v2 = v2 ./ norm(v2);
        length_a = length * norm([vec(1),vec(2)]);

        a1 = [[0.8*length_a, 0.1*length_a] * [v1;v2;] + [pt(1),pt(2)] , pt(3)+length*vec(3)]; 
        a2 = [[0.8*length_a, -0.1*length_a] * [v1;v2;] + [pt(1),pt(2)] , pt(3)+length*vec(3)];
        line([pt(1)+length*vec(1), a1(1)], [pt(2)+length*vec(2), a1(2)], [pt(3)+length*vec(3), a1(3)],'Color',c,'lineWidth',2);
        line([pt(1)+length*vec(1), a2(1)], [pt(2)+length*vec(2), a2(2)], [pt(3)+length*vec(3), a2(3)],'Color',c,'lineWidth',2);
    elseif head == 2
        % plot an arrow head on xz plane        
        v1 = [vec(1),vec(3)]; v1 = v1 ./ norm(v1);
        v2 = [vec(3),-vec(1)]; v2 = v2 ./ norm(v2);
        length_a = length * norm([vec(1),vec(3)]);

        a1 = [[0.8*length_a, 0.1*length_a] * [v1;v2;] + [pt(1),pt(3)]]; 
        a1 = [a1(1), pt(2)+length*vec(2), a1(2)];
        a2 = [[0.8*length_a, -0.1*length_a] * [v1;v2;] + [pt(1),pt(3)]];
        a2 = [a2(1), pt(2)+length*vec(2), a2(2)];
        line([pt(1)+length*vec(1), a1(1)], [pt(2)+length*vec(2), a1(2)], [pt(3)+length*vec(3), a1(3)],'Color',c,'lineWidth',2);
        line([pt(1)+length*vec(1), a2(1)], [pt(2)+length*vec(2), a2(2)], [pt(3)+length*vec(3), a2(3)],'Color',c,'lineWidth',2);  
    end
end