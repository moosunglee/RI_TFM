function Traction_vector = ...
    Get_Traction_MFISTA(h,Quiver,Young_modulus_kPa,Poisson_ratio, TV_params)

%     Quiver = Quiver_RI3; Young_modulus_kPa = 11; Poisson_ratio = 0.5; h = TFM;
%         TV_params.step = 0.005; TV_params.tv_param = 0.01; TV_params.itter_max = 100; TV_params.inner_itt = 200;

%% Shifted coordinate to prevent divergence 
    use_gpu=h.parameters.use_GPU;
% Generate a grid
    step_I = (Quiver.I(2,1,1)-Quiver.I(1,1,1));
    step_J = (Quiver.J(1,2,1)-Quiver.J(1,1,1));
    I = reshape(((1:size(Quiver.I,1)) - floor(size(Quiver.I,1)/2)-1)*step_I,[],1,1);
    J = reshape(((1:size(Quiver.I,2)) - floor(size(Quiver.I,2)/2)-1)*step_J,1,[],1);
    if size(Quiver.K,3) > 1
        step_K = (Quiver.K(1,1,2)-Quiver.K(1,1,1));
        K = reshape((0:1)*step_K,1,1,[]);
    else
        step_K = 0;
        K = step_K;
    end
    R = sqrt(I.^2+J.^2);
    
% Make a dyadic psf
    green = 2.*(1-Poisson_ratio)./(R+K) + K./R./(R+K) +...
        (2.*R.*(Poisson_ratio.*R+K)+K.^2).*I.^2./R.^3./(R+K).^2;
    green(:,:,:,1,2) = (2.*R.*(Poisson_ratio.*R+K)+K.^2).*I.*J./R.^3./(R+K).^2;
    green(:,:,:,1,3) = I.*K./R.^3 - (1-2*Poisson_ratio).*I./R./(R+K);
    green(:,:,:,2,1) = green(:,:,:,1,2);
    green(:,:,:,2,2) = 2.*(1-Poisson_ratio)./(R+K) + K./R./(R+K) +...
        (2.*R.*(Poisson_ratio.*R+K)+K.^2).*J.^2./R.^3./(R+K).^2;
    green(:,:,:,2,3) = J.*K./R.^3 - (1-2*Poisson_ratio).*J./R./(R+K);
    green(:,:,:,3,1) = I.*K./R.^3 + (1-2*Poisson_ratio).*I./R./(R+K);
    green(:,:,:,3,2) = J.*K./R.^3 + (1-2*Poisson_ratio).*J./R./(R+K);
    green(:,:,:,3,3) = K.^2./R.^3 + 2.*(1-Poisson_ratio)./R;

    green = green*(1+Poisson_ratio) / 2 / pi / (Young_modulus_kPa*1000);

    Mask = h.mk_ellipse_MS(size(green,1:2), [1 1]);
    green(isnan(green)) = 0;
    green(isinf(green)) = 0;
    green(floor(end/2)+1,floor(end/2)+1,:,:,:) = (sum(green.*Mask,1:2) / (sum(Mask(:))-1));


% Zero-pad

% FFT2
    Green = green*0;
    for j1 = 1:3
        for j2 = 1:3
            Green(:,:,:,j1,j2) = real(fftshift(fft2(ifftshift(green(:,:,:,j1,j2)))));
        end
    end

% Coordinate to draw analytic function in Fourier space

    step_Iq = 1 / step_I / size(Green,1)*2*pi;
    step_Jq = 1 / step_J / size(Green,2)*2*pi;
    Iq = reshape(((1:size(Green,1)) - floor(size(Green,1)/2)-1)*step_Iq,[],1,1);
    Jq = reshape(((1:size(Green,2)) - floor(size(Green,2)/2)-1)*step_Jq,1,[],1);
    Rq = sqrt(Iq.^2+Jq.^2);

% Green(z = 0) is analytic
    Green(:,:,1,:) = 0;
    Green(:,:,1,1,1) = (1-Poisson_ratio).*Rq.^2 + Poisson_ratio.*Jq.^2;
    Green(:,:,1,1,2) = Poisson_ratio.*Iq.*Jq;
    Green(:,:,1,2,1) = Green(:,:,1,1,2);
    Green(:,:,1,2,2) = (1-Poisson_ratio).*Rq.^2 + Poisson_ratio.*Iq.^2;
    Green(:,:,1,3,3) = (1-Poisson_ratio).*Rq.^2;
    Green(:,:,1,:) = Green(:,:,1,:) * (1 + Poisson_ratio) / pi / (Young_modulus_kPa * 1000) * 2 * pi ./ Rq.^3;
    Green(isnan(Green)) = 0;
    Green(isinf(Green)) = 0;
    Green(floor(end/2)+1,floor(end/2)+1,:,:,:) = 0;

% Inverse of Green at Z = 0 can analytically be determined.
    G_inv = Green(:,:,1,:,:)*0;
    G_inv(:,:,:,1,1) = (1-Poisson_ratio).*Rq.^2 + Poisson_ratio.*Iq.^2;
    G_inv(:,:,:,1,2) = -Poisson_ratio.*Iq.*Jq;
    G_inv(:,:,:,2,1) = G_inv(:,:,:,1,2);
    G_inv(:,:,:,2,2) = (1-Poisson_ratio).*Rq.^2 + Poisson_ratio.*Jq.^2;
    G_inv = G_inv ./ Rq ./ (1-Poisson_ratio);
    G_inv(isnan(G_inv)) = 0;
    G_inv(:,:,:,3,3) = Rq ./ (1 - Poisson_ratio);
    G_inv = G_inv*(Young_modulus_kPa*1000)/2/(1+Poisson_ratio);
    Mask = h.mk_ellipse_MS(size(green,1:2), [5 5]);
    G_inv = real(fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(G_inv))).*Mask))));


% Displacement: FFT2
    Displacement = fftshift(fft2(ifftshift(Quiver.U)));
    Displacement(:,:,:,2) = fftshift(fft2(ifftshift(Quiver.V)));
    Displacement(:,:,:,3) = fftshift(fft2(ifftshift(Quiver.W)));

% Initialize Traction vector
    Traction_vector = zeros([size(Rq,1:2),1,3],'single');
    if use_gpu
        Displacement = gpuArray(single(Displacement));
        Green = gpuArray(Green);
        G_inv = gpuArray(G_inv);
        Traction_vector = gpuArray(Traction_vector);
    end

% TV
% Initial solution of D = G * T -> ||D-GT||^2+ lambda*||T||^2
    for j1 = 1:3
        for j2 = 1:3
            Traction_vector(:,:,:,j1) = Traction_vector(:,:,:,j1) + G_inv(:,:,:,j1,j2).*Displacement(:,:,1,j2);
        end
    end
    for j1 = 1:3
        Traction_vector(:,:,:,j1) = fftshift(ifft2(ifftshift(Traction_vector(:,:,:,j1))));
    end
    Traction_vector = real(Traction_vector);
    Mask = zeros(size(Traction_vector,1:2));
    Mask(2:end-1, 2:end-1) = 1;
    Traction_vector = Traction_vector .* Mask;


% Iteration


%%

% Initialization
    err_list=[];
    dirichlet_boundary=false;
    alpha=1/TV_params.step;
    z_n=0;
    t_n=0;
    t_np=1;
    y_n=Traction_vector;
    x_n=Traction_vector;
    c_n=0;
    c_np=Inf;

    gradient_RI = Displacement*0;

    if h.parameters.verbose
        close all
        f1=figure(1);
        f2=figure(2);
        f3=figure(3);
        %{
        f4=figure(4);
        f5=figure(5);
        f6=figure(6);
        f7=figure(7);
        %}
    end
tic;
    err2 = inf;
    Dump = y_n * 0;
    stack = 0;
    for ii=1:TV_params.itter_max
        
        t_n=t_np;
        c_n=c_np;

        for j2 = 1:3
            Dump(:,:,:,j2) = fftshift(fft2(ifftshift(y_n(:,:,:,j2))));
        end

        gradient_RI(:) = 0;
        for j1 = 1:3
            for j2 = 1:3
                gradient_RI(:,:,:,j1) = gradient_RI(:,:,:,j1)+Green(:,:,:,j1,j2).*Dump(:,:,:,j2);
            end
            gradient_RI(:,:,:,j1) = fftshift(ifft2(ifftshift(gradient_RI(:,:,:,j1) - Displacement(:,:,:,j1))));
        end
        gradient_RI = real(gradient_RI);
        z_n=TV_FISTA_vector(y_n-(1/alpha)*mean(gradient_RI,3),TV_params.tv_param/alpha,...
            dirichlet_boundary,TV_params.inner_itt,use_gpu);

        % F(z_k)
        gradient_RI(:) = 0;
        for j2 = 1:3
            Dump(:,:,:,j2) = fftshift(fft2(ifftshift(z_n(:,:,:,j2))));
        end
        for j1 = 1:3
            for j2 = 1:3
                gradient_RI(:,:,:,j1) = gradient_RI(:,:,:,j1)+Green(:,:,:,j1,j2).*Dump(:,:,:,j2);
            end
            gradient_RI(:,:,:,j1) = fftshift(ifft2(ifftshift(gradient_RI(:,:,:,j1) - Displacement(:,:,:,j1))));
        end
        err = sum(abs(gradient_RI).^2,'all')/2;

        % MFISTA update
        t_np=(1+sqrt(1+4*t_n^2))/2;
        if err <= err2
            x_n=z_n;
            y_n=z_n+(t_n-1)/t_np*(z_n-x_n);
            stack = 0;
        else
            y_n=x_n+(t_n)/t_np*(z_n-x_n);
            stack = stack + 1;
        end
        Traction_vector=x_n;

        err_list = [err_list min(err,err2)];

        % F(x_(k-1))
        gradient_RI(:) = 0;
        for j2 = 1:3
            Dump(:,:,:,j2) = fftshift(fft2(ifftshift(x_n(:,:,:,j2))));
        end
        for j1 = 1:3
            for j2 = 1:3
                gradient_RI(:,:,:,j1) = gradient_RI(:,:,:,j1)+Green(:,:,:,j1,j2).*Dump(:,:,:,j2);
            end
            gradient_RI(:,:,:,j1) = fftshift(ifft2(ifftshift(gradient_RI(:,:,:,j1) - Displacement(:,:,:,j1))));
        end
        err2 = sum(abs(gradient_RI).^2,'all')/2;

        if mod(ii,20) == 0
            
            display(['itter : ' num2str(ii)]);
            if h.parameters.verbose
                set(0, 'currentfigure', f1);
                imagesc(squeeze(real(max(Traction_vector,[],3))));colorbar; axis image;title('Traction field')
                set(0, 'currentfigure', f2);
                plot(err_list);title('Cost function')
                set(0, 'currentfigure', f3);
                semilogy((err_list));title('Cost function (log)')
                %{
                set(0, 'currentfigure', f4);
                imagesc([abs(squeeze(trans_source(:,:,1,[1]))) squeeze(abs(output_field(:,:,1,[1]))) squeeze(abs(trans_source(:,:,1,1)-output_field(:,:,1,1)))]); axis image;title('Abs (predicted / experimental / delta)'),colorbar
                set(0, 'currentfigure', f5);
                imagesc([abs(squeeze(trans_source(:,:,1,[end]))) squeeze(abs(output_field(:,:,1,[scan_list(end)]))) squeeze(abs(trans_source(:,:,1,[end])-output_field(:,:,1,[scan_list(end)])))]); axis image;title('Abs (predicted / experimental / delta)'),colorbar
                set(0, 'currentfigure', f6);
                imagesc([squeeze(angle(trans_source(:,:,1,[1]))) squeeze(angle(output_field(:,:,1,[scan_list(1)]))) angle(trans_source(:,:,1,1)./output_field(:,:,1,1))]);axis image;title('Angle (predicted / experimental)'),colorbar
                set(0, 'currentfigure', f7);
                imagesc([angle(trans_source(:,:,1,end)) angle(output_field(:,:,1,scan_list(end))) angle(trans_source(:,:,1,end)./output_field(:,:,1,scan_list(end)))]);axis image;title('Angle (predicted / experimental)'),colorbar
                %}
                drawnow;
            end
        end


        if stack == 10
            break;
        end

    end
    toc;
            Traction_vector=gather(Traction_vector);

% quiver3(I(:,:,1),J(:,:,1),K(:,:,1),Traction_vector(:,:,1,1),Traction_vector(:,:,1,2),Traction_vector(:,:,1,3))
% quiver3(I(:,:,1:3),J(:,:,1:3),K(:,:,1:3),Traction_vector(:,:,1:3,1),Traction_vector(:,:,1:3,2),Traction_vector(:,:,1:3,3))
% quiver3(I(:,:,1),J(:,:,1),K(:,:,1),Quiver.U(:,:,1),Quiver.V(:,:,1),Quiver.W(:,:,1))
% quiver3(I(:,:,1:3),J(:,:,1:3),K(:,:,1:3),Quiver.U(:,:,1:3),Quiver.V(:,:,1:3),Quiver.W(:,:,1:3))
% figure,imagesc(Traction_vector(:,:,1,3)),axis image off, colormap(fireice), colorbar

% figure,imagesc(Traction_vector(:,:,1,1)),axis image off, colormap(fireice), colorbar
% figure,imagesc(Traction_vector(:,:,1,2)),axis image off, colormap(fireice), colorbar

end