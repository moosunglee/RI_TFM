function [H,JJ,II,KK] = mk_ellipse_MS(h,sizes, radii)
    if length(radii) == 1
        radii = ones(1,3) * radii;
    end
    if length(sizes) > 2
        [II,JJ,KK] = ndgrid(1:sizes(1), 1:sizes(2), 1:sizes(3));
        JJ = JJ - floor(sizes(2)/2) - 1;
        II = II - floor(sizes(1)/2) - 1;
        KK = KK - floor(sizes(3)/2) - 1;
        H = (JJ./radii(2)).^2+(II./radii(1)).^2+(KK./radii(3)).^2<=1.0;
    else
        KK = [];
        [II,JJ] = ndgrid(1:sizes(1), 1:sizes(2));
        JJ = JJ - floor(sizes(2)/2) - 1;
        II = II - floor(sizes(1)/2) - 1;
        H = (JJ./radii(2)).^2+(II./radii(1)).^2<=1.0;
    end
end