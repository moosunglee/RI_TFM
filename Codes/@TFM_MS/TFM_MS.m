classdef TFM_MS < handle
    properties (Hidden = true)
        parameters;
        utility;
    end
    methods(Static)
        function params=get_default_parameters(init_params)
            %OPTICAL PARAMETERS
%             params=BASIC_OPTICAL_PARAMETER();
            %SIMULATION PARAMETERS
            params.resolution=[1 1 1];
            params.use_GPU = true;
            % Subset spacing (pix) 
            w0 = 32;
            d0 = 6;
            params.blocksizes = [w0 w0 8]; % [iblocksize jblocksize,kblocksize]
            params.overlap = 1-d0/w0;
            params.max_shift = Inf;
            params.padding = 16;
            params.method = '3Dcentroid';
            params.N = 2;
            params.Ni = 4;
            
            if nargin==1
                params=update_struct(params,init_params);
            end
        end
    end
    methods
        function h=TFM_MS(params)
            h.parameters=params;
        end
        [imRef, imMoving] = match_size(TFM,imRef, imMoving,RI_bg);
        z_interface = Get_gel_interface(h, Im1);
        imRef = high_pass_filter(h,imRef, diameter_um);
        [C, R, cpeak] = Get_correlation_3D(h,im11, im22);
        R = peak_subpixel_positioner(h,Corr3D);
        [H,JJ,II,KK] = mk_ellipse_MS(h,sizes, radii);

        % Stitched data
        [RI_registered, shift_vector, tform] = RegisterRI_Gel(TFM,RIRef, RIMoving, z_interface_plot,RI_bg,Mask_RI_2D)
        [RI_registered, shift_vector, corr] = RegisterRI_Gel_norotation(TFM,RIRef, RIMoving, z_interface_plot,RI_bg,Mask_RI_2D)
        [RI_registered, shift_vector, tform] = RegisterFL_Gel(TFM,RIRef, RIMoving, z_interface_plot,RI_bg,Mask_RI_2D)
        Crop_mask = make_FOVmask(h, Im1, rirange,ROI_mask);
        [IZmask, JZmask] = mask_Z_mask(h, Im1, rirange)
        [RI_registered, shift_vector] = RegisterZ_Gel(TFM,RIRef, RIMoving, z_interface_plot,RI_bg,Mask_RI_2D);

        % SingleFOV
        [RI_registered, shift_vector] = RegisterRI_Gel_3D(TFM,RIRef, RIMoving, z_interface_plot,RI_bg, edge_parameter)
        [RI_registered, shift_vector, corr] = RegisterRI_Gel_Tcell_tweezers(TFM,ImRef, ImMoving, z_interface_plot,base_val,Mask_RI_2D)

        % TFM
        Quiver = Get_Quiver(TFM, imref, immoving, RI_bg, ROI_mask_final)
        Quiver = Get_Quiver_2z(TFM, ref0, moving0, ROI_mask_final)
        [Quiver, shift_vector] = QuiverShift(h,Quiver)
        Quiver = QuiverShift_v2(h,Quiver)
        [Strain, Traction_force_kPa] = Get_Strain(h,Quiver,Young_modulus_kPa,Poisson_ratio);
        Traction_vector = ...
                    Get_Traction(h,Quiver,Young_modulus_kPa,Poisson_ratio, TV_params)
        Traction_vector = ...
            Get_Traction_MFISTA(h,Quiver,Young_modulus_kPa,Poisson_ratio, TV_params)
    end
end

