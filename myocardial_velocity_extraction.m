%% Script 2. Velocity map of myocardial contraction using deformable image registration

addpath('./utils/');






%% user defined parameters
dataDim = [768, 768, 54, 158]; %row col depth time, px
voxelSize = [0.325,0.325,2];   % x, y, z direction, um
ds = [6,6,1];
dt = 5e-3;                     %time step, s
windowSize = 7;                %temporal average window size, odd integer
mask_thresholding = 300;
mask2_line_radius = 21;
load('./data2/gata1_20190718_fish4.mat');
filePath = './cmlc_gfp_gata1_dsred_control_3dpf_20190718/fish4/myocardium';
namePrefix = 'state_';
savePath = './data2/cmlc_20190718_fish4.mat';







%%
dataDim = round(dataDim./[ds,1]);
voxelSize = voxelSize.*ds;
dtt = dt*(windowSize-1);

cLine = (blood_vel.cLine+[blood_config.cropPos,0]) .* ...
    [blood_config.voxelSize(2:-1:1),blood_config.voxelSize(3)] ./ ...
    [voxelSize(2:-1:1),voxelSize(3)];  
centerPts = (blood_vel.centerPts+[blood_config.cropPos,0]) .* ...
    [blood_config.voxelSize(2:-1:1),blood_config.voxelSize(3)] ./ ...
    [voxelSize(2:-1:1),voxelSize(3)]; 

U = zeros([dataDim(1:3),dataDim(4)-windowSize+1]); V = U; W = U;
for i=ceil(windowSize/2): dataDim(4)-floor(windowSize/2)
    tic;
    
    cmlcT1 = imread3d( fullfile(filePath, [namePrefix num2str(i-floor(windowSize/2)) '.tif']) );
    cmlcT2 = imread3d( fullfile(filePath, [namePrefix num2str(i+floor(windowSize/2)) '.tif']) ); 
    cmlcT1 = imresize3(cmlcT1, dataDim(1:3));
    cmlcT2 = imresize3(cmlcT2, dataDim(1:3));
    
    [displacement,cmlcT2Reg] = imregdemons(cmlcT2,cmlcT1,[100 50 25],...
        'AccumulatedFieldSmoothing',1.3);

    mask = double( (cmlcT1 > mask_thresholding) | (cmlcT2 > mask_thresholding) );
    mask2 = distance_filter3d(mask, cLine, mask2_line_radius);

    uDemons = displacement(:,:,:,1).*mask2;  % x direction, along column
    vDemons = displacement(:,:,:,2).*mask2;  % y direction, along row
    wDemons = displacement(:,:,:,3).*mask2;  % z direction, along depth
    uDemons = uDemons * voxelSize(1) / dtt;  % convert to physical coor
    vDemons = vDemons * voxelSize(2) / dtt;
    wDemons = wDemons * voxelSize(3) / dtt;

    % visualization of vector field
    [xDemons,yDemons,zDemons] = meshgrid(1:size(mask,2),1:size(mask,1),1:size(mask,3));
    xDemons = xDemons * voxelSize(1); % x direction
    yDemons = yDemons * voxelSize(2); % y direction
    zDemons = zDemons * voxelSize(3);
    
    U(:,:,:,i-floor(windowSize/2)) = uDemons;
    V(:,:,:,i-floor(windowSize/2)) = vDemons;
    W(:,:,:,i-floor(windowSize/2)) = wDemons;
    disp(['Frame ' num2str(i) ' took ' num2str(toc) ' s']);

end

myocardium_config = struct('dataDim', dataDim, 'voxelSize', voxelSize, 'dt', dt, 'windowSize', windowSize, 'mask_thresholding', mask_thresholding,...
    'mask2_line_radius', mask2_line_radius, 'ds', ds, 'dtt', dtt);
myocardium_vel = struct('U',U,'V',V,'W',W,'X',xDemons,'Y',yDemons,'Z',zDemons,'cLine', cLine, 'centerPts', centerPts);
save(savePath, 'myocardium_config','myocardium_vel','-v7.3');

disp('Myocardial velocity map extracted.');






