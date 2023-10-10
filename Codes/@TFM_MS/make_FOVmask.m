function Crop_mask = make_FOVmask(h, Im1, rirange, Crop_mask)    
% Im1 = RIMoving.RI_stitched; rirange = [1.33 1.38];

    str_figure = 'Find XY FOV: ';
    keypressed = 0;
    if nargin == 3
        Crop_mask = Im1(:,:,1)*0;
    end
    Im1 = max(Im1,[],3);
    while keypressed ~= 1
    
    %         Im1 = edge3(Im1,"Sobel",0.3).*Im1;
    %         Im2 = edge3(Im2,"Sobel",0.3).*Im2;
    
        figure,imagesc(Im1,rirange), axis image, axis off, colormap(flip(cbrewer('seq','YlGnBu',256)))
        set(gcf,'color','w'), sgtitle([str_figure 'XY crop mask (rectangle)']),drawnow
        mask = drawrectangle;
        mask = createMask(mask);
        close
    
        Crop_mask0 = Crop_mask + mask > 0.5;
        figure, imagesc(cat(2, Im1, Crop_mask0.*Im1),rirange), axis image, colormap(flip(cbrewer('seq','YlGnBu',256)))
        set(gcf,'color','w'), title([str_figure 'Press 1 to finish the loop; otherwise retry.']),drawnow
        clc,
        keypressed = 4;
        while ~(keypressed == 1 || keypressed == 2 || keypressed == 0)
            keypressed = input([str_figure ' 1:    finish loop \n 2:    add polygon \n 3:    reset \n 0: retry:  ']);
        end
        if keypressed ~= 0
            Crop_mask = Crop_mask0;
        end
        if keypressed == 3
            Crop_mask(:)=0;
        end

        close
    end
end