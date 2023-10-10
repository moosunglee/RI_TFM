function R = peak_subpixel_positioner(h,Corr3D)
            N = 2;
           [mxv,idx] = max(Corr3D(:));
           [mi,mj,mk] = ind2sub(size(Corr3D),idx);
           if length(size(Corr3D)) > 2
                       % image near the maximum value
                    if mi-N-1>0 && mi+N+1<size(Corr3D,1) && mj-N-1>0 && mj+N+1<size(Corr3D,2) ...
                            && mk-N-1>0 && mk+N+1<size(Corr3D,3)
                    
                        im_1 = Corr3D(mi-N-1:mi+N+1,mj-N-1:mj+N+1,mk-N-1:mk+N+1);
                        % Matrix without maximum region
                        im_wm = Corr3D;
                        im_wm(mi-N-1:mi+N+1,mj-N-1:mj+N+1,mk-N-1:mk+N+1) = 0;
                        th = max(im_wm(:));
                        im_1 = im_1-th;
                        im_1(im_1<0)=0;
                        
                        meshi = mi-N-1 : mi+N+1;
                        meshj = mj-N-1 : mj+N+1;
                        meshk = mk-N-1 : mk+N+1;
                        
                        [grdi,grdj,grdk] = ndgrid(meshi,meshj,meshk); 
                        
                        ic=sum(im_1(:).*grdi(:))/sum(im_1(:));
                        jc=sum(im_1(:).*grdj(:))/sum(im_1(:));
                        kc=sum(im_1(:).*grdk(:))/sum(im_1(:));
                    else
                         ic=size(Corr3D,1)/2-1;    
                         jc=size(Corr3D,2)/2-1;    
                         kc=size(Corr3D,3)/2-1;
                    end
           else
                    if mi-N-1>0 && mi+N+1<size(Corr3D,1) && mj-N-1>0 && mj+N+1<size(Corr3D,2)
                    
                        im_1 = Corr3D(mi-N-1:mi+N+1,mj-N-1:mj+N+1);
                        im_wm = Corr3D;
                        im_wm(mi-N-1:mi+N+1,mj-N-1:mj+N+1) = 0;
                        th = max(im_wm(:));
                        % Matrix without maximum region
                        im_1(im_1<0)=0;
                        
                        meshi = mi-N-1 : mi+N+1;
                        meshj = mj-N-1 : mj+N+1;
                        
                        [grdi,grdj] = ndgrid(meshi,meshj); 
                        
                        ic=sum(im_1(:).*grdi(:))/sum(im_1(:));
                        jc=sum(im_1(:).*grdj(:))/sum(im_1(:));
                        kc=1;
                    else
                         ic=size(Corr3D,1)/2-1;    
                         jc=size(Corr3D,2)/2-1;    
                         kc=size(Corr3D,3)/2-1;
                    end
           end
            R = [ic jc kc] - floor(size(Corr3D,[1,2 3])/2) - 1;
        end