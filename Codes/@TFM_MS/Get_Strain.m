function [Strain, Traction_force_kPa] = Get_Strain(h,Quiver,Young_modulus_kPa,Poisson_ratio)
% Adapted from Hyuntae Jung


    Strain = zeros(size(Quiver.I, 1),size(Quiver.I, 2),size(Quiver.I, 3),3,3,'single');

    range0 = 1;

    ii = 1;
    while ii <= size(Quiver.I,1)
        jj = 1;
        while jj <= size(Quiver.I,2)
            kk = 1;
            while kk <= size(Quiver.I,3)
    
            % prevent 0 indices
                i1 = max(1, ii-range0*2); 
                j1 = max(1, jj-range0*2); 
                k1 = max(1, kk); 
    
            % prevent exceeding indices
                i2 = min(size(Quiver.I,1),ii+range0*2); 
                j2 = min(size(Quiver.I,2),jj+range0*2); 
                k2 = min(size(Quiver.I,3),kk+range0); 
    
            % Cut grid [um]
                I = Quiver.I(i1:i2, j1:j2, k1:k2); I = I - mean(I(:));
                J = Quiver.J(i1:i2, j1:j2, k1:k2); J = J - mean(J(:));
                K = Quiver.K(i1:i2, j1:j2, k1:k2); K = K - mean(K(:));
    
                U = Quiver.U(i1:i2, j1:j2, k1:k2);
                V = Quiver.V(i1:i2, j1:j2, k1:k2);
                W = Quiver.W(i1:i2, j1:j2, k1:k2);
                nan_list = isnan(U);
                I(nan_list) = [];
                J(nan_list) = [];
                K(nan_list) = [];
                U(nan_list) = [];
                V(nan_list) = [];
                W(nan_list) = [];
                
    
            % Get strain
                A = [I(:) J(:) K(:) ones(length(I(:)),1)];
                
                y1=(A'*A)\(A'*U(:));
                y2=(A'*A)\(A'*V(:));
                y3=(A'*A)\(A'*W(:));
    
                % Deformation gradient (https://m.blog.naver.com/richscskia/221985049233)
%                 F_matrix = [least_square_arr_u(2:4)';least_square_arr_v(2:4)';least_square_arr_w(2:4)']+eye(3);
                F_matrix = [y1(1:3) y2(1:3) y3(1:3)]' + eye(3);
                % right cauchy green tensor
                C_matrix = F_matrix.'*F_matrix;
                % Lagrangian strain tensor
                E = (C_matrix-eye(3))/2;
                Strain(ii,jj,kk,:,:) = E;
    
                kk = kk + 1;
            end
            jj = jj + 1;
        end
        ii = ii + 1;
    end
    shear_modulus = Young_modulus_kPa / (1+Poisson_ratio)/2;
    Traction_force_kPa = 2*Strain * shear_modulus;
end