function [RI_registered, shift_vector, corr] = RegisterRI_Gel_Tcell_tweezers(TFM,ImRef, ImMoving, z_interface_plot,base_val,Mask_RI_2D)

% Strategy
% 1. XY MIP translation
% 2. Crop using RIMoving
% 3. 3D correlation again


% Global registration
    if nargin == 5
        Mask_RI_2D = ImRef(:,:,1)*0+1;
    end

%% % 1. XY first
%     Mask_RI_2D_shrink = imerode(Mask_RI_2D,strel('disk',40));
    remove_edge = Mask_RI_2D*0;
    remove_edge(41:end-40,41:end-40) = 1;
%     Mask_RI_2D_shrink = remove_edge .*Mask_RI_2D_shrink;


    ImRef_XY = max(ImRef .* reshape(z_interface_plot,1,1,[]),[],3);
%     ImRef_XY = TFM.high_pass_filter(max(ImRef_XY - RI_bg,0), [0.5 0.5]);
    imMoving_XY = max(ImMoving .* reshape(z_interface_plot,1,1,[]),[],3) .* remove_edge;
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
    RI_registered = shift(ImMoving,shift_vector,base_val);
    Mask_RI_2D_shift = shift(Mask_RI_2D,[shift_vector],base_val);
%     RI_registered = imwarp(RI_registered-RI_bg,tform,'OutputView',imref2d(size(ImRef_XY)));
%     RI_registered = RI_registered + RI_bg;

%% 2. Crop using RIMoving
    [~, I1, I2, J1, J2] = cropmax(Mask_RI_2D_shift);

    ImRef_XY = ImRef(I1:I2,J1:J2,:);
    imMoving_XY = RI_registered(I1:I2,J1:J2,:);
    

% Crop
     i1 = max(1, (floor(size(ImRef_XY,1)/9)));
     j1 = max(1, (floor(size(ImRef_XY,2)/9)));
     k1 = max(1, (floor(size(ImRef_XY,3)/9)));
     i2 = min(size(ImRef_XY,1), (floor(size(ImRef_XY,1)/9*8)));
     j2 = min(size(ImRef_XY,2), (floor(size(ImRef_XY,2)/9*8)));
     k2 = min(size(ImRef_XY,3), (floor(size(ImRef_XY,3)/9*8)));
     Mask_temp = ImRef_XY*0;
     Mask_temp(i1:i2, j1:j2,k1:k2) = 1;
     ImRef_XY = ImRef_XY.*Mask_temp;
     imMoving_XY = imMoving_XY.*Mask_temp;
     ImRef_XY = ImRef_XY ./ max(ImRef_XY(:));
     imMoving_XY = imMoving_XY ./ max(imMoving_XY(:));

%%  3. Registration - 3D correlation again
    [C, shift_vector_z, corr] = TFM.Get_correlation_3D(ImRef_XY, imMoving_XY);
%     shift_vector(3) = gather(shift_vector_z(1));

% Z registration
    RI_registered = shift(RI_registered,shift_vector_z,base_val);

    shift_vector = shift_vector_z + shift_vector;

    %%
%                     im1 = max(RIRef,[],3)-RI_bg; im2 = max(RIMoving,[],3)-RI_bg; figure, imshowpair(im1, im2)
end