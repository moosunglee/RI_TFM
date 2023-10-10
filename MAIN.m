clc;clear;

% Set code path
cd0 = matlab.desktop.editor.getActiveFilename;
cd0 = cd0(1:strfind(cd0,'Codes')-2);
cds.code = fullfile(cd0, 'Codes');
addpath(genpath(cds.code));

% set dilater
dilater = strel('disk',3);
% Set result path
cds.result = fullfile(cd0);
load(fullfile(cds.result,'data.mat'))

% load TFM analyzer
resolution = [0.1834 0.1834 0.1834];
RI_bg = 1.3355;
w0 = 64; d0 = round(2/resolution(1));
params.blocksizes = [w0 w0 w0]; % [iblocksize jblocksize,kblocksize]
params.overlap = 1-d0/w0;
params.padding = 32;
params.Ni = 0;
params.N = 0;
params.use_GPU = true;
params.resolution = resolution;
params.max_shift = 450;
TFM_params=TFM_MS.get_default_parameters(params);
TFM=TFM_MS(TFM_params);

% PIV - RI
RIRef = single(RIRef) / 255*0.1145 + RI_bg;
RIMoving = single(RIMoving) / 255*0.1145 + RI_bg;
% RIRef = single(RIRef) / 40000;
% RIMoving = single(RIMoving) / 40000;
FLRef = single(FLRef);
FLMoving = single(FLMoving);
FLCell = single(FLCell);
%%

ref0 = TFM.high_pass_filter(RIRef-RI_bg, [0.5 0.5 1]).*ROI_mask_final.*ROI_mask_final2;
moving0 = TFM.high_pass_filter(RIMoving-RI_bg, [0.5 0.5 1]).*ROI_mask_final.*ROI_mask_final2;
Quiver_RI = TFM.Get_Quiver(ref0, moving0, ROI_mask_final.*ROI_mask_final2);

% PIV - FL
ref0 = FLRef.*ROI_mask_final.*ROI_mask_final2;
moving0 = FLMoving.*ROI_mask_final.*ROI_mask_final2;
Quiver_FL = TFM.Get_Quiver(ref0, moving0, ROI_mask_final.*ROI_mask_final2);


%%
    Kmax = 2;
    TV_params.step = 0.08; TV_params.tv_param = 0.025; TV_params.itter_max = 240; TV_params.inner_itt = 200;
    TFM.parameters.verbose = true;
% RI vs FL simultaneous view
% Important: U,V, W signs should be flipped
    clc,
    i1_cut = 6; 
    j1_cut = 10;
    k1_cut = 1;
    i2_cut = size(Quiver_RI.U,1)-20; 
    j2_cut = size(Quiver_RI.U,2)-6;
    k2_cut = 2;

% Trim Quiver_RI
    Quiver_RI2 = Quiver_RI;
    Quiver_RI2.U = Quiver_RI2.U(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_RI2.V = Quiver_RI2.V(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_RI2.W = Quiver_RI2.W(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_RI2.I = Quiver_RI2.I(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_RI2.J = Quiver_RI2.J(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_RI2.K = Quiver_RI2.K(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);

% Trim Quiver_FL
    Quiver_FL2 = Quiver_FL;
    Quiver_FL2.U = Quiver_FL2.U(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_FL2.V = Quiver_FL2.V(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_FL2.W = Quiver_FL2.W(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_FL2.I = Quiver_FL2.I(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_FL2.J = Quiver_FL2.J(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);
    Quiver_FL2.K = Quiver_FL2.K(i1_cut:i2_cut, j1_cut:j2_cut, k1_cut:k2_cut);

% QuiverShift
    Quiver_RI2 = TFM.QuiverShift(Quiver_RI2);
    Quiver_FL2 = TFM.QuiverShift(Quiver_FL2);

% Edge Taper
    PSF = fspecial('gaussian',min(floor(size(Quiver_RI2.J,1:2)/2)),3);
    Quiver_RI2.U = edgetaper(Quiver_RI2.U,PSF);
    Quiver_RI2.V = edgetaper(Quiver_RI2.V,PSF);
    Quiver_RI2.W = edgetaper(Quiver_RI2.W,PSF);
    Quiver_FL2.U = edgetaper(Quiver_FL2.U,PSF);
    Quiver_FL2.V = edgetaper(Quiver_FL2.V,PSF);
    Quiver_FL2.W = edgetaper(Quiver_FL2.W,PSF);

% Recalculate positions
    I1 = max(min(Quiver_RI2.I(:))-8,1);
    I2 = min(max(Quiver_RI2.I(:))+8,size(RIRef,1));
    J1 = max(min(Quiver_RI2.J(:))-8,1);
    J2 = min(max(Quiver_RI2.J(:))+8,size(RIRef,2));
    K1 = max(min(Quiver_RI2.K(:))-80,1);
    K2 = 187;

% Redesignatate positions
    Quiver_RI2.I = Quiver_RI2.I - I1 + 1;
    Quiver_RI2.J = Quiver_RI2.J - J1 + 1;
    Quiver_RI2.K = Quiver_RI2.K - K1 + 1;
    Quiver_FL2.I = Quiver_FL2.I - I1 + 1;
    Quiver_FL2.J = Quiver_FL2.J - J1 + 1;
    Quiver_FL2.K = Quiver_FL2.K - K1 + 1;


% Scale Quiver with um
    Quiver_RI2.U = -Quiver_RI2.U*resolution(1);
    Quiver_RI2.V = -Quiver_RI2.V*resolution(2);
    Quiver_RI2.W = -Quiver_RI2.W*resolution(3);
    Quiver_FL2.U = -Quiver_FL2.U*resolution(1);
    Quiver_FL2.V = -Quiver_FL2.V*resolution(2);
    Quiver_FL2.W = -Quiver_FL2.W*resolution(3);


% Prepare images to be loaded in napari
    RI = RIMoving(I1:I2, J1:J2, K1:K2);
    FLBeadAfter = FLMoving(I1:I2, J1:J2, K1:K2);
    FLBeadBefore = FLRef(I1:I2, J1:J2, K1:K2);
    FLCell = FLCell(I1:I2, J1:J2, K1:K2);


% Normalize FL images
    FLBeadAfter = max(FLBeadAfter - mean(FLBeadAfter(1:10,1:10,40:60),'all'),0);
    FLBeadAfter = round(FLBeadAfter ./ max(FLBeadAfter(:))*65535);
    FLBeadBefore = max(FLBeadBefore - mean(FLBeadBefore(1:10,1:10,40:60),'all'),0);
    FLBeadBefore = round(FLBeadBefore ./ max(FLBeadBefore(:))*65535);
    FLCell = max(FLCell - mean(FLCell(1:10,1:10,40:60),'all'),0);
    FLCell = round(FLCell ./ max(FLCell(:))*65535);

    

% Prepare traction force - RI
    Quiver_RI3 = Quiver_RI2;
    Quiver_RI3.I = Quiver_RI3.I*resolution(1);
    Quiver_RI3.J = Quiver_RI3.J*resolution(2);
    Quiver_RI3.K = Quiver_RI3.K*resolution(3);
    Quiver_RI3.U = Quiver_RI3.U(:,:,1:Kmax);
    Quiver_RI3.V = Quiver_RI3.V(:,:,1:Kmax);
    Quiver_RI3.W = Quiver_RI3.W(:,:,1:Kmax);
    Quiver_RI3.I = Quiver_RI3.I(:,:,1:Kmax);
    Quiver_RI3.J = Quiver_RI3.J(:,:,1:Kmax);
    Quiver_RI3.K = Quiver_RI3.K(:,:,1:Kmax);
    if  isempty(dir('Traction_RI.mat*'))
        Traction_RI = TFM.Get_Traction(Quiver_RI3,11,0.5, TV_params);
        save('Traction_RI.mat', 'Traction_RI');
    else
        load('Traction_RI.mat')
    end
        
% Prepare traction force - FL
    Quiver_FL3 = Quiver_FL2;
    Quiver_FL3.I = Quiver_FL3.I*resolution(1);
    Quiver_FL3.J = Quiver_FL3.J*resolution(2);
    Quiver_FL3.K = Quiver_FL3.K*resolution(3);
    Quiver_FL3.U = Quiver_FL3.U(:,:,1:Kmax);
    Quiver_FL3.V = Quiver_FL3.V(:,:,1:Kmax);
    Quiver_FL3.W = Quiver_FL3.W(:,:,1:Kmax);
    Quiver_FL3.I = Quiver_FL3.I(:,:,1:Kmax);
    Quiver_FL3.J = Quiver_FL3.J(:,:,1:Kmax);
    Quiver_FL3.K = Quiver_FL3.K(:,:,1:Kmax);

    if  isempty(dir('Traction_FL.mat*'))
        Traction_FL = TFM.Get_Traction(Quiver_FL3,11,0.5, TV_params);

        save('Traction_FL.mat', 'Traction_FL');
    else
        load('Traction_FL.mat')
    end

% Prepare intepolation grid
    interpolater_II = 1:size(RI,1);
    interpolater_JJ = 1:size(RI,2);
    [II, JJ] = ndgrid(interpolater_II, interpolater_JJ);

% Get traction - RI 
    V = Traction_RI(:,:,1,1); 
    U = Traction_RI(:,:,1,2);
    W = Traction_RI(:,:,1,3);
    I = round(Quiver_RI3.I(:,:,1)/resolution(1),3)-8; 
    J = round(Quiver_RI3.J(:,:,1)/resolution(1),3)-8;
    K = round(Quiver_RI3.K(:,:,1)/resolution(1))+5-24;
    I_RI = I(2:end-1,2:end-1);
    J_RI = J(2:end-1,2:end-1);
    K_RI = K(2:end-1,2:end-1);
    U_RI = U(2:end-1,2:end-1);
    V_RI = V(2:end-1,2:end-1);
    W_RI = W(2:end-1,2:end-1);

% Get traction - FL
    V = Traction_FL(:,:,1,1); 
    U = Traction_FL(:,:,1,2);
    W = Traction_FL(:,:,1,3);
    I = round(Quiver_FL3.I(:,:,1)/resolution(1),3)-8; 
    J = round(Quiver_FL3.J(:,:,1)/resolution(1),3)-8;
    K = round(Quiver_FL3.K(:,:,1)/resolution(1))+5;
    I_FL = I(2:end-1,2:end-1);
    J_FL = J(2:end-1,2:end-1);
    K_FL = K(2:end-1,2:end-1);
    U_FL = U(2:end-1,2:end-1);
    V_FL = V(2:end-1,2:end-1);
    W_FL = W(2:end-1,2:end-1);


%%
% Draw Fig. 01b
cmap_green = gray(255).*[0 1 0];
ri_range = [1.333 1.37];
fl_range = [0 40000];

img = repmat(FLBeadAfter, [1 1 1 3])/65535*10;
crop0 = 7500; fl_mag = 3;
img(:,:,:,2) = img(:,:,:,2) + max(FLCell-crop0,0)/65535*fl_mag;

img = min(img, 65535);

figure,imagesc(squeeze(img(:,:,91,:))),axis image off, title('Z = 92 (2 um)'), hold on
V = Traction_RI(:,:,1,1); 
U = Traction_RI(:,:,1,2);
W = Traction_RI(:,:,1,3);

cmap = cbrewer('div', 'RdBu',255);
cmap(cmap<0) = 0;
vectorMagnitude = hypot(U(:),V(:)); 
vecMagNorm = min(vectorMagnitude / 1000,1);
Colormagnitude = (max(min(W(:) / 500,1),-1)+1)/2;
vecColorIdx = round(Colormagnitude * (size(cmap,1)-1)) + 1; 
scale_factor = 0.05;


for i = 1:numel(U)
    if vecMagNorm(i) > 0.1
        quiver(J(i),I(i),U(i)*scale_factor,V(i)*scale_factor,...
            'Autoscale','off','LineWidth',(10)*vecMagNorm(i),'ShowArrowHead','on','MaxHeadSize', (10)*vecMagNorm(i),...
            'color', cmap(vecColorIdx(i),:) )
        hold on
    end
end
axis image off, set(gcf,'color','w')
hold off
drawnow;


figure,imagesc(squeeze(img(:,:,91,:))),axis image off, title('Z = 92 (2 um)'), hold on
V = Traction_FL(:,:,1,1); 
U = Traction_FL(:,:,1,2);
W = Traction_FL(:,:,1,3);
    % im0=imagesc(max(img(9:end-8,9:end-8,1:K(1)-6),[],3));axis off,axis image,axis off, colormap(gray)

cmap = cbrewer('div', 'RdBu',255);
cmap(cmap<0) = 0;
vectorMagnitude = hypot(U(:),V(:)); 
vecMagNorm = min(vectorMagnitude / 1000,1);
Colormagnitude = (max(min(W(:) / 500,1),-1)+1)/2;
vecColorIdx = round(Colormagnitude * (size(cmap,1)-1)) + 1; 


for i = 1:numel(U)
    if vecMagNorm(i) > 0.1
        quiver(J(i),I(i),U(i)*scale_factor,V(i)*scale_factor,...
            'Autoscale','off','LineWidth',(10)*vecMagNorm(i),'ShowArrowHead','on','MaxHeadSize', (10)*vecMagNorm(i),...
            'color', cmap(vecColorIdx(i),:) )
        hold on
    end
end
axis image off, set(gcf,'color','w')
hold off
drawnow;


%%

fl_range = [0 40000];

img = repmat(FLBeadAfter, [1 1 1 3])/65535*10;
crop0 = 7500; fl_mag = 3;
img(:,:,:,2) = img(:,:,:,2) + max(FLCell-crop0,0)/65535*fl_mag;

img = min(img, 65535);

% figure,imagesc(squeeze(img(:,:,91,:))),axis image off, title('Z = 92 (2 um)'), hold on
figure,imagesc(RI(:,:,79),ri_range),axis image off, title('Z = 92 (2 um)'), colormap((gray)), hold on
V = Quiver_RI3.U(:,:,1); 
U = Quiver_RI3.V(:,:,1);
W = Quiver_RI3.W(:,:,1);
    % im0=imagesc(max(img(9:end-8,9:end-8,1:K(1)-6),[],3));axis off,axis image,axis off, colormap(gray)

cmap = cbrewer('div', 'RdBu',255);
cmap(cmap<0) = 0;
vectorMagnitude = hypot(U(:),V(:)); 
vecMagNorm = min(vectorMagnitude / max(vectorMagnitude),1);
Colormagnitude = (max(min(W(:) / max(vectorMagnitude),1),-1)+1)/2;
vecColorIdx = round(Colormagnitude * (size(cmap,1)-1)) + 1; 
scale_factor = 40;


for i = 1:numel(U)
    if vecMagNorm(i) > 0.1
        quiver(J(i),I(i),U(i)*scale_factor,V(i)*scale_factor,...
            'Autoscale','off','LineWidth',(5)*vecMagNorm(i),'ShowArrowHead','on','MaxHeadSize', (5)*vecMagNorm(i),...
            'color', 'r' )
        hold on
    end
end
axis image off, set(gcf,'color','w')
hold off
drawnow;



figure,imagesc(squeeze(img(:,:,91,:))),axis image off, title('Z = 92 (2 um)'), hold on
V = Quiver_FL3.U(:,:,1); 
U = Quiver_FL3.V(:,:,1);
W = Quiver_FL3.W(:,:,1);
    % im0=imagesc(max(img(9:end-8,9:end-8,1:K(1)-6),[],3));axis off,axis image,axis off, colormap(gray)

cmap = cbrewer('div', 'RdBu',255);
cmap(cmap<0) = 0;
vectorMagnitude = hypot(U(:),V(:)); 
vecMagNorm = min(vectorMagnitude / max(vectorMagnitude),1);
Colormagnitude = (max(min(W(:) / max(vectorMagnitude),1),-1)+1)/2;
vecColorIdx = round(Colormagnitude * (size(cmap,1)-1)) + 1; 


for i = 1:numel(U)
    if vecMagNorm(i) > 0.1
        quiver(J(i),I(i),U(i)*scale_factor,V(i)*scale_factor,...
            'Autoscale','off','LineWidth',(5)*vecMagNorm(i),'ShowArrowHead','on','MaxHeadSize', (5)*vecMagNorm(i),...
            'color', 'c' )
        hold on
    end
end
axis image off, set(gcf,'color','w')
hold off
drawnow;



sqrt(mean(abs(Quiver_RI3.U(:,:,1:2)-Quiver_FL3.U(:,:,1:2)).^2,'all'))*resolution(3)*1000
sqrt(mean(abs(Quiver_RI3.V(:,:,1:2)-Quiver_FL3.V(:,:,1:2)).^2,'all'))*resolution(3)*1000
sqrt(mean(abs(Quiver_RI3.W(:,:,1:2)-Quiver_FL3.W(:,:,1:2)).^2,'all'))*resolution(3)*1000


%%
close all
     figure, im0=imagesc(permute(squeeze(max(RI(9:end-8,9:end-8,25:end-30,:),[],1)),[2,1,3]),ri_range);axis off,axis image,axis off, colormap gray
%         figure(404),  im3=imagesc(squeeze(max(RI_Moving(9:end-8,9:end-8,1:end,:),[],1))', [1.34 1.39]);axis off,axis image, colormap((gray))
        hold on
        J = squeeze(Quiver_RI3.J(1,:,1:2)/resolution(1)); K = squeeze(Quiver_RI3.K(1,:,1:2)/resolution(1))-24+5;
        V = absmax(Quiver_RI3.U(:,:,1:2),1); 
        U = absmax(Quiver_RI3.V(:,:,1:2),1);
        W = absmax(Quiver_RI3.W(:,:,1:2),1);
        U = U(2:end-1);
        V = V(2:end-1);
        W = W(2:end-1);

        % Assign colors based on magnitude of vectors
         vectorMagnitude = hypot(U(:),W(:)); 
        vecMagNorm = min(vectorMagnitude / max(vectorMagnitude),1);
        Colormagnitude = (max(min(W(:) / 500,1),-1)+1)/2;
        vecColorIdx = round(Colormagnitude * (size(cmap,1)-1)) + 1; 

        for i = 1:numel(V)
            if vecMagNorm(i) > 0.01
                quiver(J(i),K(i),U(i)*scale_factor,W(i)*scale_factor,...
                    'Autoscale','off','LineWidth',(5)*vecMagNorm(i),'ShowArrowHead','on','MaxHeadSize', (5)*vecMagNorm(i),...
                    'color', 'c' )
                hold on
            end
        end
        axis image off, set(gcf,'color','w')
        hold off



img = repmat(max(FLBeadAfter-500,0), [1 1 1 3])/65535*3;
crop0 = 6000; fl_mag = 2.1;
img(:,:,:,2) = img(:,:,:,2) +max(FLCell-crop0,0)/65535*fl_mag;
img = min(img, 1);

img = min(img, 65535);
     % z traction
     figure, im0=imagesc(permute(squeeze(max(img(9:end-8,9:end-8,25:end-30,:),[],1)),[2,1,3]));axis off,axis image,axis off, 
%         figure(404),  im3=imagesc(squeeze(max(RI_Moving(9:end-8,9:end-8,1:end,:),[],1))', [1.34 1.39]);axis off,axis image, colormap((gray))
        hold on
        J = squeeze(Quiver_FL3.J(1,:,1:2)/resolution(1)); K = squeeze(Quiver_FL3.K(1,:,1:2)/resolution(1))-24+5;
        V = absmax(Quiver_FL3.U(:,:,1:2),1); 
        U = absmax(Quiver_FL3.V(:,:,1:2),1);
        W = absmax(Quiver_FL3.W(:,:,1:2),1);
        U = U(2:end-1);
        V = V(2:end-1);
        W = W(2:end-1);

        % Assign colors based on magnitude of vectors
         vectorMagnitude = hypot(U(:),W(:)); 
        vecMagNorm = min(vectorMagnitude / max(vectorMagnitude),1);
        Colormagnitude = (max(min(W(:) / 500,1),-1)+1)/2;
        vecColorIdx = round(Colormagnitude * (size(cmap,1)-1)) + 1; 

        for i = 1:numel(V)
            if vecMagNorm(i) > 0.01
                quiver(J(i),K(i),U(i)*scale_factor,W(i)*scale_factor,...
                    'Autoscale','off','LineWidth',(5)*vecMagNorm(i),'ShowArrowHead','on','MaxHeadSize', (5)*vecMagNorm(i),...
                    'color', 'm' )
                hold on
            end
        end
        axis image off, set(gcf,'color','w')
        hold off


