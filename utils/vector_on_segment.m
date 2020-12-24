function [pointer, magnitude, location] = vector_on_segment(U,V,W,X,Y,Z,M,M_segment,segment_ind)


    u = U(M_segment==segment_ind);v = V(M_segment==segment_ind);w = W(M_segment==segment_ind);
    x = X(M_segment==segment_ind);y = Y(M_segment==segment_ind);z = Z(M_segment==segment_ind);
    s = M(M_segment==segment_ind);
    pointer = (s' * [u,v,w]) ./ sum(s);
    pointer = pointer ./ norm(pointer);
    magnitude = mean(s);
    location = [x, y, z];
    
end