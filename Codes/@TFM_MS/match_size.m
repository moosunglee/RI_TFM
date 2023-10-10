function [imRef, imMoving] = match_size(TFM,imRef, imMoving,RI_bg)
    
        % Set 3D volume size of PIV
                size_3D = max(size(imRef), size(imMoving));
            
        % Padding RIRef.RI_multi
            if norm(size_3D - size(imRef)) ~=0
                imRef = My_paddzero_MS(imRef, size_3D, RI_bg);
            end

        % Padding RIMoving.RI_multi
            if norm(size_3D - size(imMoving)) ~=0
                imMoving = My_paddzero_MS(imMoving, size_3D, RI_bg);
            end   
    
end