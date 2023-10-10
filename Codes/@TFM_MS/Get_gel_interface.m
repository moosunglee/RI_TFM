function z_interface = Get_gel_interface(h, Im1)

        str_figure = 'Get z interface';
        figure, imagesc(squeeze(max(Im1,[],1))'), axis image, axis off, colormap(flip(cbrewer('seq','YlGnBu',256)))
        set(gcf,'color','w'), sgtitle([str_figure ' (XZ crop mask)']),drawnow
        mask = drawpolygon;
        mask = createMask(mask);

        close

        [~, i1,i2,j1,j2] = cropmax(mask);

        z_interface = i1;

end
