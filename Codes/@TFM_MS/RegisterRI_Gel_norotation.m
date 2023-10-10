function [RI_registered, shift_vector, corr] = RegisterRI_Gel_norotation(TFM,RIRef, RIMoving, z_interface_plot,RI_bg,Mask_RI_2D)

% Global registration
    if nargin == 5
        Mask_RI_2D = RIRef(:,:,1)*0+1;
    end
    Mask_RI_2D_shrink = imerode(Mask_RI_2D,strel('disk',40));
    remove_edge = Mask_RI_2D_shrink*0;
    remove_edge(41:end-40,41:end-40) = 1;
    Mask_RI_2D_shrink = remove_edge .*Mask_RI_2D_shrink;

% 1. XY first


    ImRef_XY = max(RIRef .* reshape(z_interface_plot,1,1,[]),[],3) .* Mask_RI_2D_shrink;
%     ImRef_XY = TFM.high_pass_filter(max(ImRef_XY - RI_bg,0), [0.5 0.5]);
    imMoving_XY = max(RIMoving .* reshape(z_interface_plot,1,1,[]),[],3) .* remove_edge;
%     imMoving_XY = TFM.high_pass_filter(max(imMoving_XY - RI_bg,0), [0.5 0.5]);

% Crop
     i1 = max(1, (floor(size(ImRef_XY,1)/9)));
     j1 = max(1, (floor(size(ImRef_XY,2)/9)));
     i2 = min(size(ImRef_XY,1), (floor(size(ImRef_XY,1)/9*8)));
     j2 = min(size(ImRef_XY,2), (floor(size(ImRef_XY,2)/9*8)));
     Mask_temp = ImRef_XY*0;
     Mask_temp(i1:i2, j1:j2) = 1;
     ImRef_XY = ImRef_XY.*Mask_temp;
     imMoving_XY = imMoving_XY.*Mask_temp;
     ImRef_XY = ImRef_XY ./ max(ImRef_XY(:));
     imMoving_XY = imMoving_XY ./ max(imMoving_XY(:));

% Registration - (1) translation (2) rotation
    [C, shift_vector, corr] = TFM.Get_correlation_3D(ImRef_XY, imMoving_XY);
    corr
    shift_vector = gather(shift_vector);
%     imMoving_XY = shift(imMoving_XY,shift_vector,0);
%     [optimizer, metric] = imregconfig('multimodal');
%     optimizer.InitialRadius = 0.009;
%     optimizer.Epsilon = 1.5e-2;
%     optimizer.GrowthFactor = 1.01;
%     optimizer.MaximumIterations = 300;
%     remove_edge = remove_edge * 0;
%     remove_edge(101:end-100,101:end-100) = 1;
%     tform = imregtform(imMoving_XY, ImRef_XY, 'rigid', optimizer, metric);
%     movingRegistered = imwarp(imMoving_XY,tform,'OutputView',imref2d(size(ImRef_XY)));
%                     figure, imshowpair(ImRef_XY/0.01,movingRegistered/0.01)
% % %                     figure, imagesc((C)), axis image off, drawnow

% XY registration
    RI_registered = shift(RIMoving,[shift_vector 0],RI_bg);
%     RI_registered = imwarp(RI_registered-RI_bg,tform,'OutputView',imref2d(size(ImRef_XY)));
%     RI_registered = RI_registered + RI_bg;

% 2. XZ next
    ImRef_XZ = squeeze(max(RIRef.* Mask_RI_2D_shrink .* reshape(z_interface_plot,1,1,[]),[],1))' ;
%     ImRef_XZ = TFM.high_pass_filter(max(ImRef_XZ - RI_bg,0), [1 0.5]);
    imMoving_XZ = squeeze(max(RI_registered.* Mask_RI_2D_shrink,[],1))';
%     imMoving_XZ = TFM.high_pass_filter(max(imMoving_XZ - RI_bg,0), [1 0.5]);
%                     figure, imshowpair(ImRef_XZ,imMoving_XZ)

% Crop
     i1 = max(1, (floor(size(ImRef_XZ,1)/9)));
     j1 = max(1, (floor(size(ImRef_XZ,2)/9)));
     i2 = min(size(ImRef_XZ,1), (floor(size(ImRef_XZ,1)/9*8)));
     j2 = min(size(ImRef_XZ,2), (floor(size(ImRef_XZ,2)/9*8)));
     Mask_temp = ImRef_XZ*0;
     Mask_temp(i1:i2, j1:j2) = 1;
     ImRef_XZ = ImRef_XZ.*Mask_temp;
     imMoving_XZ = imMoving_XZ.*Mask_temp;
     ImRef_XZ = ImRef_XZ ./ max(ImRef_XZ(:));
     imMoving_XZ = imMoving_XZ ./ max(imMoving_XZ(:));

% Registration - (1) translation (2) rotation
    [C, shift_vector_z, corr] = TFM.Get_correlation_3D(ImRef_XZ, imMoving_XZ);
    shift_vector(3) = gather(shift_vector_z(1));

% Z registration
    RI_registered = shift(RI_registered,[0 0 shift_vector(3)],RI_bg);

    %%
%                     im1 = max(RIRef,[],3)-RI_bg; im2 = max(RIMoving,[],3)-RI_bg; figure, imshowpair(im1, im2)
end