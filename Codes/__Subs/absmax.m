function img2 = absmax(img,dim)

    img2 = max(abs(img),[],dim);

    idxx = find(abs(img) ~= img2);
    img(idxx) = -inf;

    img2(:) = max(img,[],dim);

end