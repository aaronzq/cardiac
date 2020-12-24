function strain_rate = compute_strain_rate_from_line(coor1, coor2, U, V, W, M, X, Y, Z, mask_percentile, mask_thickness, voxel_size)
% coor1: coordinate of line, x, y, z
% coor2: coordinate of line, x, y, z
% mask_percentile: 
% mask_thickness: px
    debug = 0;

    dcoor = coor2-coor1;
    dir = find(abs(dcoor)==max(abs(dcoor)));
    sample_line = zeros(max(abs(dcoor))+1,3);
    sample_linePhy = sample_line;    
    if(dcoor(dir)) < 0
        temp = coor1;
        coor1 = coor2;
        coor2 = temp;
        dcoor = coor2-coor1;
    end   
    for i = coor1(dir):coor2(dir) 
        x = coor1(1) + (i-coor1(dir))/(coor2(dir)-coor1(dir))*(coor2(1)-coor1(1));
        y = coor1(2) + (i-coor1(dir))/(coor2(dir)-coor1(dir))*(coor2(2)-coor1(2));
        z = coor1(3) + (i-coor1(dir))/(coor2(dir)-coor1(dir))*(coor2(3)-coor1(3));
        sample_line(i-coor1(dir)+1,:) = round( [y,x,z]./voxel_size );
        sample_linePhy(i-coor1(dir)+1,:) = [y,x,z];
    end
    
    mask0 = ones(size(X));
    mask1 = distance_filter3d(mask0,sample_line(1:round(max(dcoor)*mask_percentile),:),mask_thickness);
    mask2 = distance_filter3d(mask0,sample_line(end-round(max(dcoor)*mask_percentile):end,:),mask_thickness);


    sample_vels = zeros(2,size(M,4));
    sample_distance = zeros(1, size(M,4));
    for t=1:size(M,4)

        u = U(:,:,:,t);
        v = V(:,:,:,t);  % y, row
        w = W(:,:,:,t);
        m = M(:,:,:,t);  % magnitude

        u1 = u(find(mask1==1));v1 = v(find(mask1==1));w1 = w(find(mask1==1));
        x1 = X(find(mask1==1));y1 = Y(find(mask1==1));z1 = Z(find(mask1==1));
        s1 = m(find(mask1==1));
        pointer1 = (s1' * [u1,v1,w1]) ./ sum(s1);
        pointer1Pos = ((s1>0)' * [x1,y1,z1]) ./ sum(s1>0);
        pointer1_on_line = dcoor / norm(dcoor) * pointer1';

        u2 = u(find(mask2==1));v2 = v(find(mask2==1));w2 = w(find(mask2==1));
        x2 = X(find(mask2==1));y2 = Y(find(mask2==1));z2 = Z(find(mask2==1));
        s2 = m(find(mask2==1));
        pointer2 = (s2' * [u2,v2,w2]) ./ sum(s2);
        pointer2Pos = ((s2>0)' * [x2,y2,z2]) ./ sum(s2>0);
        pointer2_on_line = dcoor / norm(dcoor) * pointer2';

        sample_vels(1,t) = pointer1_on_line;
        sample_vels(2,t) = pointer2_on_line;
        sample_distance(t) = norm(pointer2Pos-pointer1Pos);
        if debug
            if t==1
                figure;
                quiver3(Y,X,Z,v,u,w,4,'Color',[0.1,0.5,0.3],'LineWidth', 0.7);hold on;
                quiver3(Y,X,Z,v.*mask1,u.*mask1,w.*mask1,2,'Color',[0.9,0.1,0.1],'LineWidth', 0.9);
                quiver3(Y,X,Z,v.*mask2,u.*mask2,w.*mask2,2,'Color',[0.9,0.1,0.1],'LineWidth', 0.9);
                scatter3(sample_linePhy(:,1), sample_linePhy(:,2), sample_linePhy(:,3),25,[0.1 0.1 0.1],'filled');hold off;
                axis equal
                view([90 90]);  
                xlim([0,250]); ylim([0,250]); zlim([0,200]); 
                set(gca, 'Projection', 'perspective'); 
                set(gca, 'GridColor', 'k');
                set(gca, 'lineWidth', 1);
                set(gca, 'GridAlpha', 0.3);
                set(gcf,'Color','w');
                grid on; 
            end
        end
   
    end

    strain_rate = (sample_vels(2,:) - sample_vels(1,:))/max(sample_distance);



end