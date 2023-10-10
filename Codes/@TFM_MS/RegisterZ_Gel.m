function [RI_registered, shift_vector] = RegisterZ_Gel(TFM,RIRef, RIMoving, z_interface_plot,RI_bg,Mask_RI_2D)

% Global registration
    if nargin == 5
        Mask_RI_2D = RIRef(:,:,1)*0+1;
    end
    Mask_RI_2D_shrink = imerode(Mask_RI_2D,strel('disk',40));
    remove_edge = Mask_RI_2D_shrink*0;
    remove_edge(41:end-40,41:end-40) = 1;
    Mask_RI_2D_shrink = remove_edge .*Mask_RI_2D_shrink;


% 2. XZ next
    ImRef_XZ = squeeze(max(RIRef.* Mask_RI_2D_shrink .* reshape(z_interface_plot,1,1,[]),[],1))' ;
%     ImRef_XZ = TFM.high_pass_filter(max(ImRef_XZ - RI_bg,0), [1 0.5]);
    imMoving_XZ = squeeze(max(RIMoving.* Mask_RI_2D_shrink,[],1))';
%     imMoving_XZ = TFM.high_pass_filter(max(imMoving_XZ - RI_bg,0), [1 0.5]);
%                     figure, imshowpair(ImRef_XZ,imMoving_XZ)

% Crop
    shift_vector = [0 0 0];
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
    RI_registered = shift(RIMoving,[0 0 shift_vector(3)],RI_bg);

    %%
%                     im1 = max(RIRef,[],3)-RI_bg; im2 = max(RIMoving,[],3)-RI_bg; figure, imshowpair(im1, im2)
end