%% Script 1. Velocity map of bloof flow from Imaris automatic tracking results
%  Coordinate system: [x,y,z] corresponds to [col,row,depth] in the image matrix
clear all;
addpath('./utils/');





%% User defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%% dataset parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDim = [76, 39, 54, 158];  %row col depth time, px
cropPos = [13, 30];           %crop position in original image, row col, px
voxelSize = [2, 2, 2];        %x, y, z direction, um
dt = 5e-3;                    %time step, s
debug = true;                %debug modes prints all intermediate figures
%%%%%%%%%%%%%%%%%%%%%%%%%%% center line parameter %%%%%%%%%%%%%%%%%%%%%%%%%
cLineResolution = 1000;       %number of pts on center line
cLineAdj1 = 0;                %offset at the beginning of the line, rad
cLineAdj2 = 0;                %offset at the end of the line, rad
invert = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%% file path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracks = importTrackImaris('./data2/gata1_20190718_fish4_Position.csv');
saveName = './data2/gata1_20190718_fish4.mat';





%% Extract vector maps
nTracks = size(tracks,2);
U = zeros(dataDim); V = U; W = U;
xgg=0;
for i=1:nTracks
    
    track = tracks{i};
    trackImpute = zeros(max(track(:,1))-min(track(:,1))+1, size(track,2) );
    trackImputeFrame = min(track(:,1)) : max(track(:,1));
    % Data imputation due to the potential skip steps (if maxmimum tracking step larger than 1 in Imaris)
    for t = 1:size(trackImputeFrame,2)
        ind = find(track(:,1)==trackImputeFrame(t));
        if isempty(ind)
            indd = find(track(:,1)==trackImputeFrame(t-1));
            trackImpute(t,:) = 0.5*(track(indd,:) + track(indd+1,:));
        else
            trackImpute(t,:) = track(ind,:);
        end
    end
    trackImputeVel = [];
    trackImputeVel(:,1) = trackImpute(:,1); % frame
    trackImputeVel(:,2:4) = round(trackImpute(:,2:end) ./ voxelSize); % positions in pixel: x, y, z
    temp = trackImpute(2:end,2:end) - trackImpute(1:end-1,2:end); 
    trackImputeVel(:,5:7) = ([temp;0,0,0]+[0,0,0;temp]) ./ [1,1,1; repmat(2,[size(temp,1)-1,3]); 1,1,1]; % velocity in um: u, v, w
    trackImputeVel(:,5:7) = trackImputeVel(:,5:7) ./ dt; % velocity in um/s: u(x-direction), v(y-direction), w(z-direction); 
    for v = 1:size(trackImputeVel,1)
        xgg=xgg+1;
        U(trackImputeVel(v,3), trackImputeVel(v,2), trackImputeVel(v,4), trackImputeVel(v,1)) = trackImputeVel(v,5); %row col depth time
        V(trackImputeVel(v,3), trackImputeVel(v,2), trackImputeVel(v,4), trackImputeVel(v,1)) = trackImputeVel(v,6);
        W(trackImputeVel(v,3), trackImputeVel(v,2), trackImputeVel(v,4), trackImputeVel(v,1)) = trackImputeVel(v,7);
    end
end
% trackImputeVel: row(px), col(px), depth(px), frame(px), vel_x(um/s), vel_y(um/s), vel_z(um/s)
% U: vel_x(um/s)
% V: vel_y(um/s)
% W: vel_z(um/s)
[xRBC,yRBC,zRBC] = meshgrid(1:size(U(:,:,:,t),2),1:size(U(:,:,:,t),1),1:size(U(:,:,:,t),3));
xRBC = xRBC*voxelSize(1);    yRBC = yRBC*voxelSize(2);    zRBC = zRBC*voxelSize(3);
M = zeros(size(U));
for t = 1:size(U,4)
    M(:,:,:,t) = sqrt(U(:,:,:,t).^2 + V(:,:,:,t).^2 + W(:,:,:,t).^2); 
end
clear i ind t temp track trackImpute trackImputeFrame trackImputeVel v xgg

%% PCA analysis on the orientation of blood cells distribution
UPts = (sum(U~=0,4)>0 | sum(V~=0,4)>0 | sum(W~=0,4)>0);
[UPtsCorr1,UPtsCorr2,UPtsCorr3] = ind2sub(size(UPts),find(UPts>0));  %UPtsCorr1: row indice;  UPtsCorr2: col indice; UPtsCorr3: depth indice
[V_,S, D_] = pca([UPtsCorr1,UPtsCorr2,UPtsCorr3]);
v1 = V_(:,1); v2 = V_(:,2); v3 = V_(:,3);
pt0 = mean([UPtsCorr1,UPtsCorr2,UPtsCorr3]);
if debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    zoom = 0.16;
    pt1 = [ pt0(1)+v1(1)*D_(1)*zoom, pt0(2)+v1(2)*D_(1)*zoom, pt0(3)+v1(3)*D_(1)*zoom ];
    pt2 = [ pt0(1)+v2(1)*D_(2)*zoom, pt0(2)+v2(2)*D_(2)*zoom, pt0(3)+v2(3)*D_(2)*zoom ];
    pt3 = [ pt0(1)+v3(1)*D_(3)*zoom, pt0(2)+v3(2)*D_(3)*zoom, pt0(3)+v3(3)*D_(3)*zoom ];
    scatter3(UPtsCorr1,UPtsCorr2,UPtsCorr3);hold on;
    line([pt0(1),pt1(1)],[pt0(2),pt1(2)],[pt0(3),pt1(3)],'Color','red','LineWidth',2);
    line([pt0(1),pt2(1)],[pt0(2),pt2(2)],[pt0(3),pt2(3)],'Color','blue','LineWidth',2);
    line([pt0(1),pt3(1)],[pt0(2),pt3(2)],[pt0(3),pt3(3)],'Color','magenta','LineWidth',2);
    axis equal; hold off;
    title('PCA analysis on blood cells distribution.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Cells sorting on the projection plane
P = [v2,v1]; % projection matrix
UPtsCorr = [UPtsCorr1,UPtsCorr2,UPtsCorr3];
UPtsCorr_P = UPtsCorr * P;
centerPts_P = round([max(UPtsCorr_P(:,1))+1, 0.5*(max(UPtsCorr_P(:,2)) + min(UPtsCorr_P(:,2)))]);
pt0_P = pt0 * P;
centerPts = pt0 + (centerPts_P-pt0_P) * P';
UVecs_P = [UPtsCorr_P(:,1)-centerPts_P(1), UPtsCorr_P(:,2)-centerPts_P(2)];
UVecs_P(:,3) = atan( UVecs_P(:,2) ./ UVecs_P(:,1) );
[s,ind] = sort(UVecs_P(:,3));
UVecs_P_sort = UVecs_P(ind,:);
UPtsCorr_P_sort = UPtsCorr_P(ind,:);
UPtsCorr_sort = UPtsCorr(ind,:);

if debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c = linspace(1,10,size(UPtsCorr_sort,1));
    figure;
    scatter3(UPtsCorr_sort(:,1),UPtsCorr_sort(:,2),UPtsCorr_sort(:,3),25,c', 'filled'); axis equal; hold on;
    scatter3(centerPts(1),centerPts(2),centerPts(3), 25, 'k', 'filled'); hold off;
    colorbar
    title('Sort the pts on original space')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

clear zoom pt0 pt0_P pt1 pt2 pt3 pline1 pline2 pline3 V_ S D_
%% Fitting a center line through blood cells
theta = UVecs_P_sort(:,3);
y1 = UPtsCorr_sort(:,1); y2 = UPtsCorr_sort(:,2); y3 = UPtsCorr_sort(:,3);
TT = [theta.^4, theta.^3, theta.^2, theta, ones(size(theta))];
p1 = TT\y1; p2 = TT\y2; p3 = TT\y3;
theta_fit = linspace(min(theta)+cLineAdj1,max(theta)+cLineAdj2,cLineResolution);
theta_fit = theta_fit';
if invert
    theta_fit = theta_fit(end:-1:1);
end
TT_fit = [theta_fit.^4, theta_fit.^3, theta_fit.^2, theta_fit, ones(size(theta_fit))];
UPts_fit = [TT_fit * p1, TT_fit * p2, TT_fit * p3];

cLine = [UPts_fit(:,1), UPts_fit(:,2), UPts_fit(:,3)];
[cLineV1, cLineV2, cLineV3, cLineV11] = axis_on_line(cLine, centerPts);

if debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    scatter3(UPtsCorr_sort(:,1),UPtsCorr_sort(:,2),UPtsCorr_sort(:,3),25,c', 'filled'); axis equal; hold on;
    line(cLine(:,1),cLine(:,2),cLine(:,3), 'Color', 'k', 'lineStyle', '-', 'LineWidth', 3); 
    scatter3(centerPts(1),centerPts(2),centerPts(3), 25, 'k', 'filled'); 
    cc1 = [repmat(centerPts(1), 1, cLineResolution); cLine(:,1)'];
    cc2 = [repmat(centerPts(2), 1, cLineResolution); cLine(:,2)'];
    cc3 = [repmat(centerPts(3), 1, cLineResolution); cLine(:,3)'];
    line(cc1(:), cc2(:), cc3(:), 'Color', 'r');
    for iii=1:10:cLineResolution
    plot_vec([cLine(iii,1),cLine(iii,2),cLine(iii,3)],cLineV1(iii,:),10,'g',0);
    plot_vec([cLine(iii,1),cLine(iii,2),cLine(iii,3)],cLineV2(iii,:),10,'m',0);
    plot_vec([cLine(iii,1),cLine(iii,2),cLine(iii,3)],cLineV3(iii,:),10,'c',0);
    end
    hold off;
    colorbar;
    title('Center line and reference coor system');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
clear centerPts_P cvx* p1 p2 p3 s theta TT UPts UPts_fit UPtsCorr_P UPtsCorr_P_sort UPtsCorr1 UPtsCorr2 UPtsCorr3 UPtsCorr UVecs* v1 v2 v3 y1 y2 y3 P ind c
clear cc1 cc2 cc3 cLineAdj1 cLineAdj2 iii invert temp theta_fit TT_fit 
%% Save data
blood_config = struct('dataDim', dataDim, 'voxelSize', voxelSize, 'dt', dt, 'cLineResolution', cLineResolution,'cropPos',cropPos);
blood_vel = struct('U',U,'V',V,'W',W,'M',M,'X',xRBC,'Y',yRBC,'Z',zRBC,'UPtsCorr_sort', UPtsCorr_sort,...
    'centerPts',centerPts,'cLine',cLine,'cLineV1',cLineV1,'cLineV2',cLineV2,'cLineV3',cLineV3);
save(saveName, 'blood_config','blood_vel','-v7.3');
disp('Blood flow velocity map extracted.');








