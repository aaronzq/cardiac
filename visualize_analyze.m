%% Script 3. Data analysis and visualization
clear all;
addpath('./utils');
load('./data2/cmlc_20190718_fish4.mat')
load('./data2/gata1_20190718_fish4.mat')

%% Plot blood cell velocity map
savePath0 = './temp2/blood_velocity';
if ~exist(savePath0, 'dir')
    mkdir(savePath0);
end

clear vels1d meanVels1D
blood_vel.cLine_phy = blood_vel.cLine * diag([blood_config.voxelSize(2:-1:1),blood_config.voxelSize(3)]);
for t=ceil(myocardium_config.windowSize/2):size(blood_vel.U,4)-floor(myocardium_config.windowSize/2)  % in order to sycn with cmlc vec map, which starts from 4 to (end-3)
    tic;
    tt = t-floor(myocardium_config.windowSize/2);
    vels1d(:,tt) = vector_map_project_onto_line1d(blood_vel.U(:,:,:,t),blood_vel.V(:,:,:,t),blood_vel.W(:,:,:,t),blood_vel.cLine_phy,blood_config.voxelSize);
    if sum(vels1d(:,tt)~=0,1) > 0
        meanVels1D(tt) = sum(vels1d(:,tt),1) ./ sum(vels1d(:,tt)~=0,1);
    else
        meanVels1D(tt) = 0;
    end    
    figure;
    quiver3(blood_vel.Y,blood_vel.X,blood_vel.Z,blood_vel.V(:,:,:,t),blood_vel.U(:,:,:,t),blood_vel.W(:,:,:,t),40,'Color',[0.99,0.3,0.3],...
        'LineWidth', 2); hold on;
    scatter3(blood_vel.cLine_phy(:,1),blood_vel.cLine_phy(:,2),blood_vel.cLine_phy(:,3),10 ,'k', 'filled'); hold off; 
    axis equal 
    view([110 76])
    % set(gca, 'Ydir', 'reverse')
    xlim([0,170]); ylim([0,120]); zlim([0,200]);
    set(gca, 'Projection', 'perspective');
    set(gca, 'GridColor', 'k');
    set(gca, 'lineWidth', 1);
    set(gca, 'GridAlpha', 0.3);
    set(gcf,'Color','w');
    grid on;
    title(['Time ' num2str(tt*blood_config.dt*1000) 'ms']);      
%     print(gcf, '-dpng', '-r150', fullfile(savePath0, ['RBC_' num2str(tt) '.png']));
    saveas(gcf, fullfile(savePath0, ['RBC_' num2str(tt) '.png']));
    close(gcf);

    disp(['Frame ' num2str(tt) ' took ' num2str(toc) ' s']);
end

for t = 1:size(vels1d,2)
    flowDistance(t) = blood_config.dt*sum( meanVels1D(1:t), 'all' );
end

figure;
subplot(2,1,1); plot((1:size(vels1d,2))*blood_config.dt, meanVels1D,'k','LineWidth',1);
title('Flow velocity (um/s)');
subplot(2,1,2); plot((1:size(vels1d,2))*blood_config.dt, flowDistance,'k','LineWidth',1);
title('Flow distance (um)');
xlabel('Time (s)');
set(gcf,'Color','w');

%% Plot myocardial velocity map
ds = [3, 3, 3]; % spatial downsample rate for vector map visualization
cLineRange = 1:1000;
savePath1 = 'temp2/myocardial_velocity_xy';
savePath2 = 'temp2/myocardial_velocity_yz';


if ~exist(savePath1, 'dir')
    mkdir(savePath1);
end
if ~exist(savePath2, 'dir')
    mkdir(savePath2);
end
[M,globalM] = magnitude_of_vels(myocardium_vel.U, myocardium_vel.V, myocardium_vel.W);
myocardium_vel.cLine_phy = myocardium_vel.cLine * diag([myocardium_config.voxelSize(2:-1:1),myocardium_config.voxelSize(3)]);
myocardium_vel.centerPts_phy = myocardium_vel.centerPts * diag([myocardium_config.voxelSize(2:-1:1),myocardium_config.voxelSize(3)]);

contractionDirection1 = zeros([size(M,4),1]);
contractionDirection2 = contractionDirection1; 
contractionDirection3 = contractionDirection1;
contractionDirection4 = contractionDirection1;
contractionVelocity1 = zeros([size(M,4),1]);
contractionVelocity2 = contractionVelocity1;
contractionVelocity3 = contractionVelocity1;
contractionVelocity4 = contractionVelocity1;
for t = 1 : size(M,4)
    tic;   
    M_segment = divide_segments(M(:,:,:,t), myocardium_vel.cLine_phy(cLineRange,:), myocardium_vel.centerPts_phy, myocardium_config.voxelSize);
    [vector1, magnitude1, location1] = vector_on_segment(myocardium_vel.U(:,:,:,t), myocardium_vel.V(:,:,:,t), myocardium_vel.W(:,:,:,t),...
        myocardium_vel.X, myocardium_vel.Y, myocardium_vel.Z, M(:,:,:,t), M_segment, 1);
    [vector2, magnitude2, location2] = vector_on_segment(myocardium_vel.U(:,:,:,t), myocardium_vel.V(:,:,:,t), myocardium_vel.W(:,:,:,t),...
        myocardium_vel.X, myocardium_vel.Y, myocardium_vel.Z, M(:,:,:,t), M_segment, 2);
    [vector3, magnitude3, location3] = vector_on_segment(myocardium_vel.U(:,:,:,t), myocardium_vel.V(:,:,:,t), myocardium_vel.W(:,:,:,t),...
        myocardium_vel.X, myocardium_vel.Y, myocardium_vel.Z, M(:,:,:,t), M_segment, 3);
    [vector4, magnitude4, location4] = vector_on_segment(myocardium_vel.U(:,:,:,t), myocardium_vel.V(:,:,:,t), myocardium_vel.W(:,:,:,t),...
        myocardium_vel.X, myocardium_vel.Y, myocardium_vel.Z, M(:,:,:,t), M_segment, 4);

    contractionDirection1(t) = compute_angle_from_vector(vector1(1), -vector1(2), 0);
    contractionDirection2(t) = compute_angle_from_vector(vector2(1), -vector2(2), 0);
    contractionDirection3(t) = compute_angle_from_vector(vector3(1), -vector3(2), 0);
    contractionDirection4(t) = compute_angle_from_vector(vector4(1), -vector4(2), 0);
    contractionVelocity1(t) = magnitude1;
    contractionVelocity2(t) = magnitude2;
    contractionVelocity3(t) = magnitude3;
    contractionVelocity4(t) = magnitude4;    
    
    if t == 41       
        figure;
        scatter3(location1(:,2),location1(:,1),location1(:,3),25,[0.8 0.1 0.1],'filled'); hold on;
        scatter3(location2(:,2),location2(:,1),location2(:,3),25,[0.1 0.6 0.1],'filled');
        scatter3(location3(:,2),location3(:,1),location3(:,3),25,[0.1 0.1 0.8],'filled');
        scatter3(location4(:,2),location4(:,1),location4(:,3),25,[0.9 0.5 0.3], 'filled'); hold off;
        axis equal;
        set(gca, 'Projection', 'orthographic');
        set(gca, 'GridColor', 'k');
        set(gca, 'lineWidth', 1);
        set(gca, 'GridAlpha', 0.3);
        set(gcf,'Color','w');
        xlim([0,250]); ylim([0,250]); zlim([0,200]); 
        axis off       
    end
    
    U_ds = myocardium_vel.U(1:ds(1):end, 1:ds(2):end, 1:ds(3):end, t);
    V_ds = myocardium_vel.V(1:ds(1):end, 1:ds(2):end, 1:ds(3):end, t);
    W_ds = myocardium_vel.W(1:ds(1):end, 1:ds(2):end, 1:ds(3):end, t);
    X_ds = myocardium_vel.X(1:ds(1):end, 1:ds(2):end, 1:ds(3):end);
    Y_ds = myocardium_vel.Y(1:ds(1):end, 1:ds(2):end, 1:ds(3):end);
    Z_ds = myocardium_vel.Z(1:ds(1):end, 1:ds(2):end, 1:ds(3):end);

    figure;
    quiver3(Y_ds,X_ds,Z_ds,V_ds,U_ds,W_ds,4,'Color',[0.1,0.5,0.3],'LineWidth', 0.7);hold on;
    plot_vec(mean([location1(:,2) location1(:,1) location1(:,3)]), [vector1(2) vector1(1) vector1(3)], 30, 'r', 1);
    plot_vec(mean([location2(:,2) location2(:,1) location2(:,3)]), [vector2(2) vector2(1) vector2(3)], 30, 'r', 1);
    plot_vec(mean([location3(:,2) location3(:,1) location3(:,3)]), [vector3(2) vector3(1) vector3(3)], 30, 'r', 1);
    plot_vec(mean([location4(:,2) location4(:,1) location4(:,3)]), [vector4(2) vector4(1) vector4(3)], 30, 'r', 1);
    hold off;
    axis equal
    view([90 90]);  
    xlim([0,250]); ylim([0,250]); zlim([0,200]); 
    set(gca, 'Projection', 'perspective'); 
    set(gca, 'GridColor', 'k');
    set(gca, 'lineWidth', 1);
    set(gca, 'GridAlpha', 0.3);
    set(gcf,'Color','w');
    grid on;          
    title(['Time ' num2str(t*myocardium_config.dt*1000) 'ms']); 
    saveas(gcf, fullfile(savePath1, ['demons_xy_' num2str(t) '.png']));
    close(gcf);
    figure;
    quiver3(Y_ds,X_ds,Z_ds,V_ds,U_ds,W_ds,4,'Color',[0.1,0.5,0.3],'LineWidth', 0.7);hold on;
    plot_vec(mean([location1(:,2) location1(:,1) location1(:,3)]), [vector1(2) vector1(1) vector1(3)], 30, 'r', 2);
    plot_vec(mean([location2(:,2) location2(:,1) location2(:,3)]), [vector2(2) vector2(1) vector2(3)], 30, 'r', 2);
    plot_vec(mean([location3(:,2) location3(:,1) location3(:,3)]), [vector3(2) vector3(1) vector3(3)], 30, 'r', 2);
    plot_vec(mean([location4(:,2) location4(:,1) location4(:,3)]), [vector4(2) vector4(1) vector4(3)], 30, 'r', 2);
    hold off;   
    axis equal
    view([-180,0]); 
    xlim([0,250]); ylim([0,250]); zlim([0,200]); 
    set(gca, 'Projection', 'perspective'); 
    set(gca, 'GridColor', 'k');
    set(gca, 'lineWidth', 1);
    set(gca, 'GridAlpha', 0.3);
    set(gcf,'Color','w');
    grid on;   
    saveas(gcf, fullfile(savePath2, ['demons_yz_' num2str(t) '.png']));
    close(gcf);   

    
    disp(['Frame ' num2str(t) ' took ' num2str(toc) ' s']);   
end
figure;
plot( (1:size(M,4))*myocardium_config.dt, angle_unwrap(contractionDirection1),...
    'Color', [0.6 0.2 0.2], 'LineWidth', 2); hold on
plot( (1:size(M,4))*myocardium_config.dt, angle_unwrap(contractionDirection2),...
    'Color', [0.2 0.6 0.2], 'LineWidth', 2);
plot( (1:size(M,4))*myocardium_config.dt, angle_unwrap(contractionDirection3),...
    'Color', [0.2 0.2 0.6], 'LineWidth', 2);
plot( (1:size(M,4))*myocardium_config.dt, angle_unwrap(contractionDirection4),...
    'Color', [0.9 0.5 0.3], 'LineWidth', 2); 
hold off;
title('Segmental contractile directions (unwrapped, degree)');
xlabel('Time (s)');
set(gcf,'Color','w');


figure;
plot( (1:size(M,4))*myocardium_config.dt, contractionVelocity1,...
    'Color', [0.6 0.2 0.2],'LineWidth', 2); hold on
plot( (1:size(M,4))*myocardium_config.dt, contractionVelocity2,...
    'Color', [0.2 0.6 0.2], 'LineWidth', 2);
plot( (1:size(M,4))*myocardium_config.dt, contractionVelocity3,...
    'Color', [0.2 0.2 0.6], 'LineWidth', 2);
plot( (1:size(M,4))*myocardium_config.dt, contractionVelocity4,...
    'Color', [0.9 0.5 0.3], 'LineWidth', 2); 
plot( (1:size(M,4))*myocardium_config.dt, globalM,...
    'k--', 'LineWidth', 2);
hold off;
title('Segmental contractile magnitudes (um/s)');
xlabel('Time (s)');
set(gcf,'Color','w');

%% Compute strain rate from velocity map 

% define a line in 3D
coor1 = [28,58,80];
coor2 = [169,58,80];

% parameters used in sampling
mask_percentile = 0.5;
mask_thickness = 2;

% sample the vectors and compute strain rate
strain_rate = compute_strain_rate_from_line(coor1, coor2, myocardium_vel.U, myocardium_vel.V,...
    myocardium_vel.W, M, myocardium_vel.X, myocardium_vel.Y, myocardium_vel.Z, mask_percentile,...
    mask_thickness, myocardium_config.voxelSize);

figure;
plot((1:length(strain_rate))*myocardium_config.dt, strain_rate,'k-', 'LineWidth', 2);
title('Strain rate');
xlabel('Time (s)');
set(gcf,'Color','w');     
    
    
