function [RI_registered, shift_vector] = RegisterRI_Gel_3D(TFM,RIRef, RIMoving, z_interface_plot,RI_bg, edge_parameter)
% Global registration
    if nargin == 5
        edge_parameter = 20;
    end
    Mask = zeros(size(RIRef),'single');
    Mask(edge_parameter+1:end-edge_parameter,edge_parameter+1:end-edge_parameter) = 1;
    RIRef = RIRef .* Mask  .* reshape(z_interface_plot,1,1,[]);
    RIMoving = RIMoving .* Mask  .* reshape(z_interface_plot,1,1,[]);

% Registration - (1) translation (2) rotation
    [C, shift_vector, corr] = TFM.Get_correlation_3D(RIRef, RIMoving);
    shift_vector = gather(shift_vector);
    corr
% XY registration
    RI_registered = shift(RIMoving,shift_vector,RI_bg);
    RI_registered = RI_registered + RI_bg;

    %%
%                     im1 = max(RIRef,[],3)-RI_bg; im2 = max(RIMoving,[],3)-RI_bg; figure, imshowpair(im1, im2)
%                     im1 = squeeze(max(RIRef,[],1))'-RI_bg; im2 = squeeze(max(RIMoving,[],1))'-RI_bg; figure, imshowpair(im1, im2)
% figure,orthosliceViewer(C)

end