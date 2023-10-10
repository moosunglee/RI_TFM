function Mask = compress_mask(mask)

    Mask = [];

    Mask{1} = max(mask,[],1);
    Mask{2} = max(mask,[],2);
    Mask{3} = max(mask,[],3);

end