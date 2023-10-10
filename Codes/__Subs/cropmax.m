function    [RI_B2, i1,i2,j1,j2,k1,k2] = cropmax(RI_B)
    [xind,yind,zind]=ind2sub(size(RI_B),find(RI_B>0.1));
    i1=min(xind);j1=min(yind);k1=min(zind);
   i2=max(xind);j2=max(yind);k2=max(zind);
%     x1=min(x1,y1);x2=max(x2,y2);
%     y1=x1;y2=x2;
    
    RI_B2=RI_B(i1:i2,j1:j2,k1:k2);
end