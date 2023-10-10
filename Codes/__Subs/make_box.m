function mask0 = make_box(mask0)

%         CC = bwconncomp(mask0, 26);
%         numPixels = cellfun(@numel,CC.PixelIdxList);
%         [biggest,idx0] = max(numPixels);
%         for jj = 1:length(CC.PixelIdxList)
%             if jj ~= idx0
%                 mask0(CC.PixelIdxList{jj}) = 0;
%             end
%         end


        [~, i1, i2, j1, j2,k1,k2] = cropmax(mask0);
        mask0(:)=0;
        mask0(i1:i2, j1:j2, k1:k2) = 1;
end