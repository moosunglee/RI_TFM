function Quiver = Get_Quiver_2z(TFM, ref0, moving0, ROI_mask_final)
% ROI_mask_final = ROI_mask .* RIRef.Mask_RI_2D .* reshape(z_interface_plot, 1, 1, []);

if nargin == 3
    ROI_mask_final = ones(size(ref0),'single');
end

[ROI_mask_final, i1, i2, j1, j2, k1, k2] = cropmax(ROI_mask_final);

if (i1 - TFM.parameters.blocksizes(1)/2>=1) && (j1 - TFM.parameters.blocksizes(2)/2>=1) &&...
        (k1 - TFM.parameters.blocksizes(3)/2>=1) && (i2 + TFM.parameters.blocksizes(1)/2<=size(ROI_mask_final,1)) &&...
        (j2 + TFM.parameters.blocksizes(2)/2<=size(ROI_mask_final,2)) && (k2 + TFM.parameters.blocksizes(3)/2<=size(ROI_mask_final,3))

    ref0 = ref0(i1 - TFM.parameters.blocksizes(1)/2:i2 + TFM.parameters.blocksizes(1)/2,...
        j1 - TFM.parameters.blocksizes(2)/2:j2 + TFM.parameters.blocksizes(2)/2,...
        k1 - TFM.parameters.blocksizes(3)/2:k2 + TFM.parameters.blocksizes(3)/2);
    ref0(k1 - TFM.parameters.blocksizes(3)/2:k2:k1) = 0;
    moving0 = moving0(i1 - TFM.parameters.blocksizes(1)/2:i2 + TFM.parameters.blocksizes(1)/2,...
        j1 - TFM.parameters.blocksizes(2)/2:j2 + TFM.parameters.blocksizes(2)/2,...
        k1 - TFM.parameters.blocksizes(3)/2:k2 + TFM.parameters.blocksizes(3)/2);
    moving0(k1 - TFM.parameters.blocksizes(3)/2:k2:k1) = 0;
else
    ref0 = ref0(i1:i2, j1:j2, k1:k2);
    moving0 = moving0(i1:i2, j1:j2, k1:k2);
    ref0 = padarray(ref0, TFM.parameters.blocksizes/2,0);
    moving0 = padarray(moving0, TFM.parameters.blocksizes/2,0);
end

% figure,orthosliceViewer(ROI_mask_final),pause
ROI_mask_final = padarray(ROI_mask_final, TFM.parameters.blocksizes/2,0);


%         imshowpair(squeeze(max(ref0,[],1)), squeeze(max(moving0,[],1)))
        if size(ref0,3) < TFM.parameters.blocksizes(3)
            TFM.parameters.blocksizes(3) = size(ref0,3);
        end

    % GPU settings
        if TFM.parameters.use_GPU
            ref0 = gpuArray(single(ref0));
            moving0 = gpuArray(single(moving0));
            ROI_mask_final = gpuArray(ROI_mask_final);
        end
        
        if any(size(ref0) ~= size(moving0))
            error('Image dimensions do not match together.')
        end
    
    % 3. Define tile increments
        incs = round(TFM.parameters.blocksizes* (1-TFM.parameters.overlap)); %[inci, incj, inck] - define increment
        if sum((incs < 1) + (incs > TFM.parameters.blocksizes)) > 0 
            error('Wrong Overlap in Correlation Algorithm')
        end
        
    % 4. Initialize Quiver
        Quiver = struct;
        sizes = size(moving0);
        if size(ref0,3) > 1
            Qsz = floor((sizes) ./ incs); % If it meets the ceiling, an exceptional coordinate will be used
        else
            Qsz = floor((sizes(1:2)) ./ incs(1:2));
            Qsz(3) = 1;
        end
        if Qsz(3) > 1
            Qsz(3) = 2;
        end


        Quiver.U = zeros(Qsz,'single');
        Quiver.V = zeros(Qsz,'single');
        Quiver.W = zeros(Qsz,'single');
        if length(sizes(3)) == 2
            Quiver.K = zeros(Qsz,'single')+inf;
            Quiver.W = zeros(Qsz,'single')+inf;
        end
        Quiver.pkh = zeros(Qsz,'single');


%
    f = waitbar(0,'Get Quiver');
    ki = 1; kj = 1; kk = 1;
    t1 = clock;
    for ki = 1:Qsz(1)
%         clc,
%         disp([num2str(ki) ' / ' num2str(Qsz(1))]);
        t2 = etime(clock,t1);
        waitbar(ki / Qsz(1),f,['PIV is in progress, ' num2str(t2) ' seconds...']);
        for kj = 1:Qsz(2)
            for kk = 1:Qsz(3)

            % initialize iterative process
                niterations=0; DI=inf; DJ=inf; DK=inf; 
                % Criterion: when sub-pixel registration meets certain lower bound
                while niterations<=TFM.parameters.N && abs(Quiver.U(ki, kj, kk)-DI) > 0.02 &&...
                        abs(Quiver.V(ki, kj, kk)-DJ) > 0.02 && abs(Quiver.W(ki, kj, kk)-DK) > 0.02
                % Allocate the displacement on the new matrix
                    niterations = niterations + 1;
                    DI = Quiver.U(ki, kj, kk);
                    DJ = Quiver.V(ki, kj, kk);
                    DK = Quiver.W(ki, kj, kk);

                % Define subwindow index
                    idxs_start = ([ki kj kk]-1).*incs+1;
                    idxs_end = idxs_start + TFM.parameters.blocksizes-1;
                    idxs_end = min(idxs_end, size(ref0,1:3));

                % Get ref0 and imMoving

                % Get shifted window coordinate
                    subpixel_shift = [0 0 0];
                    if ~(DI == 0 && DJ == 0 && DK == 0) %skip first iteration(DX=0,Dy=0 for first iteration)
                        idxs_start2 = idxs_start + round([DI,DJ,DK]);
                        idxs_start2 = max(min(idxs_start2, sizes), [1,1,1]);
                        idxs_end2 = idxs_end + round([DI,DJ,DK]);
                        idxs_end2 = max(min(idxs_end2, sizes), [1,1,1]);
                        subpixel_shift = [DI,DJ,DK] - round([DI,DJ,DK]);
                        idxs_end = idxs_start + (idxs_end2-idxs_start2);
                    else
                        idxs_start2 = idxs_start;
                        idxs_end2 = idxs_end;
                    end

                    if Qsz(3) > 1
                        imRef = ref0(idxs_start(1):idxs_end(1), idxs_start(2):idxs_end(2), idxs_start(3):idxs_end(3));
                        check_zero = ROI_mask_final(idxs_start(1):idxs_end(1), idxs_start(2):idxs_end(2), idxs_start(3):idxs_end(3));
                    else
                        imRef = ref0(idxs_start(1):idxs_end(1), idxs_start(2):idxs_end(2), :);
                        check_zero = ROI_mask_final(idxs_start(1):idxs_end(1), idxs_start(2):idxs_end(2), :);
                    end
                % Get ref0 and imMoving
                    if Qsz(3) > 1
                        imMoving = moving0(idxs_start2(1):idxs_end2(1), idxs_start2(2):idxs_end2(2), idxs_start2(3):idxs_end2(3));
                    else
                        imMoving = moving0(idxs_start2(1):idxs_end2(1), idxs_start2(2):idxs_end2(2),:);
                    end
                    if norm(subpixel_shift) ~=0
                        imMoving = shift(imMoving,subpixel_shift,0);
                    end

                    if sum(check_zero(:)) > 0
                    % zeropadding subwindow
                        imRef = padarray(imRef, ones(1,length(sizes))*TFM.parameters.padding,0);
                        imMoving = padarray(imMoving, ones(1,length(sizes))*TFM.parameters.padding,0);
    
                    % PIV
                        [~, dshift, Quiver.pkh(ki,kj,kk)] = TFM.Get_correlation_3D(imRef, imMoving);

%                         clc,dshift
    
                    % End of computation
                        Quiver.U(ki, kj, kk) = Quiver.U(ki, kj, kk) + dshift(1);
                        Quiver.V(ki, kj, kk) = Quiver.V(ki, kj, kk) + dshift(2);
                        Quiver.W(ki, kj, kk) = Quiver.W(ki, kj, kk) + dshift(3);
                    end
%% End of loop
                end
            end
        end
    end
close(f)
    %%

        Quiver.U(isnan(Quiver.U)) = 0;
        Quiver.V(isnan(Quiver.V)) = 0;
        Quiver.W(isnan(Quiver.W)) = 0;

    % Get index
        II = ((1:Qsz(1))-1)*incs(1)+1+TFM.parameters.blocksizes(1)/2;
        JJ = ((1:Qsz(2))-1)*incs(2)+1+TFM.parameters.blocksizes(2)/2;
        KK = ((1:Qsz(3))-1)*incs(3)+1+TFM.parameters.blocksizes(3)/2;
    
        II(end) = min(II(end), ((Qsz(1)-1)*incs(1)+sizes(1))/2);
        JJ(end) = min(JJ(end), ((Qsz(2)-1)*incs(2)+sizes(2))/2);
        KK(end) = min(KK(end), ((Qsz(3)-1)*incs(3)+sizes(3))/2);
        II = II - TFM.parameters.blocksizes(1)/2;
        JJ = JJ - TFM.parameters.blocksizes(2)/2;
        KK = KK - TFM.parameters.blocksizes(2)/2;
    
        [Quiver.I, Quiver.J, Quiver.K] = ndgrid(II,JJ,KK);
        Quiver.J = Quiver.J + j1 -1;
        Quiver.I = Quiver.I + i1 - 1;
        Quiver.K = Quiver.K + k1 -1;


% 
% 
% scale_factor = 5;
% figure,
% im0=imagesc(max(immoving(:,:,1:end),[],3), [1.336 1.38]);axis off,axis image,axis off, colormap(flip(gray))
%             hold on
%             I = mean(Quiver.I,3); J = mean(Quiver.J,3);
%             V = absmax(Quiver.V(:,:,4:6),3); U = absmax(Quiver.U(:,:,4:6),3);
%             quiver(J,I,...
%                         -(V)*scale_factor, -(U)*scale_factor,'magenta','Autoscale','off','ShowArrowHead','on'), axis image, set(gcf,'color','w'), axis off
% 
% 
%             figure,
%             im3=imagesc(squeeze(max(immoving,[],1))', [1.336 1.38]);axis off,axis image, colormap(flip(gray))
%             hold on
%             J = squeeze(Quiver.J(1,:,4:end-5)); K = squeeze(Quiver.K(1,:,4:end-5));
%             V = squeeze(absmax(Quiver.V(floor(end/2)-3:floor(end/2)+3,:,4:end-5),1)); 
%             W =  squeeze(absmax(Quiver.W(floor(end/2)-3:floor(end/2)+3,:,4:end-5),1)); 
%             quiver(J(:),K(:),...
%                         (V(:))*scale_factor, (W(:))*scale_factor,'magenta','Autoscale','off','ShowArrowHead','on'), axis image, set(gcf,'color','w'), axis off
%   
% 
