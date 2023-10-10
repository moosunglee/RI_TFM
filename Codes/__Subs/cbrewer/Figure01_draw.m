clc;clear;close all
gpuDevice(1)
reset(gpuDevice)
%% ODT_solver_MS
cdcode = 'C:\Users\Moosung Lee\Desktop\Moosung_CBS'; 
% cdcode = mfilename('fullpath'); idxstr = strfind(cdcode, '\'); cdcode = cdcode(1:idxstr(end)-1);
cd(cdcode)
addpath(genpath(cdcode))
% cddata = 'D:\_CODES\_Hardwares\_Labview\00000Hardware rearrangement v3\Meadowlark SLM\TEMP\Tetraspeck_100nm_sp\set000002'; 
% generalvars = {'GPU_on', 'cdsp'}; parameters = {1, cddata};
% h = ODT_solver_MS(generalvars,parameters);
h = ODT_solver_MS_v2_clean;
% load_data_ODT(h,h.cdbg_sp,h.cdsp);
% ODT_spec_to_CBS(h) % - Under construction


%% Initialize parameters
%%%%% params_optics
params_optics = struct;
params_optics.lambda = 0.532;
params_optics.n_m = 1.336;
params_optics.dx = params_optics.lambda / params_optics.n_m / 8;
params_optics.use_GPU = true;
params_optics.use_cuda = true;
params_optics.return_3D = true;
params_optics.return_transmission = true;
params_optics.return_reflection = true; %ensure the 1rst order reflection is fully simulated (might take longer)
params_optics.verbose = false;
params_optics.ZP = round(26.*params_optics.lambda/params_optics.n_m/params_optics.dx.*[0.26 0.26 0.52]); % [# of voxels]
% scan_NA =  linspace(0.6,1,3);
params_optics.obj_NA = 1.2;
params_optics.scan_NA = params_optics.n_m;

%%%%% Sample
Sample = struct;
Sample.name = 'Sphere';
Sample.radius = 1.25; % Sample.radius = 2.5; 1.5; % [um; bead Sample.radius]
Sample.n_s = 1.4; % Sample.n_s = 1.5; % Sample.n_s = 1.4;%815; % RI of PMMA @ 1064 nm % Sample.n_s = 1.4496; % RI of fused silica @ 1064 nm % Sample.n_s = 1.5718; % RI of PS @clc 1050 nm
% Sample.n_s = 0.1215+params_optics.n_m;
Sample.antialiasing_n = 1;
h.params_optics = params_optics; h.Sample = Sample; clear Sample params_optics



%%%%% incident plane wave parameters
%{
pols = [1; 0]; k0s = [0; 0]; 
%}
% %{
scan_NA = 1.08;
thetas = linspace(0, 2*pi, 21);
thetas(end) = [];
pols = 1;
k0s = [0; 0];
%}


generate_sample(h); % Make a RImap
h.Sample.RImap = circshift(h.Sample.RImap, [0 0 ceil(h.Sample.radius / h.params_optics.dx)])+ circshift(h.Sample.RImap, -ceil([0 0 h.Sample.radius / h.params_optics.dx]))-h.params_optics.n_m;
% h.Sample.RImap = h.Sample.RImap * 0 + h.params_optics.n_m;
%%
% Make_isotropic_tensor(h) % To test dyadic computation % orthosliceViewer(gather(real(RImap(:,:,:,1,1))))
generate_source(h,pols,k0s); % Prepare incident plane waves   
clear pols k0s

% h.Sample.RImap = circshift(h.Sample.RImap, [5 10 3]);

%% forward
clc;
MSE = []; Peak_SE = [];
j1_list = 1;
for j1 = j1_list
h.source.xtol = 1e-8;
h.source.padding_source = 0;
h.source.sink_abs = 0.3;

h.source.iterations_number = -1;
h.params_optics.verbose = true;
h.params_optics.verbose_display_period = 100;
% h.params_optics.Method = 'Convergent Born';
% h.source.boundary_mode = 'ACC';
% h.source.sink_thickness = 10; % [lambda/n_m
h.params_optics.Method = 'Convergent Born - analytic';
h.source.boundary_mode = 'ACC';
h.source.sink_thickness = [1.5 1.5 3]; % [lambda/n_m

% 
% h.params_optics.Method = 'Convergent Born';
% h.source.boundary_mode = 'Vanilla';
% h.source.sink_thickness = 26;

[h, Field3D, Field2D_trans, Field2D_ref] = forward(h);
h.Results.Field2D_trans = Field2D_trans;
h.Results.Field3D = Field3D;
h.Results.Field2D_ref = Field2D_ref;
clear Field2D_trans Field3D Field2D_ref
end
%%
crange = [-3 3];
subplot(141),imagesc(real(h.Results.Field3D(:,:,end-6)),crange),axis image, axis off
subplot(142),imagesc(real(h.Results.Field3D(:,:,end-4)),crange),axis image, axis off
subplot(143),imagesc(real(h.Results.Field3D(:,:,end-2)),crange),axis image, axis off
subplot(144),imagesc(real(h.Results.Field3D(:,:,end-0)),crange),axis image, axis off
% subplot(141),imagesc(real(h.Results.Field3D(:,:,end-6))),axis image, axis off
% subplot(142),imagesc(real(h.Results.Field3D(:,:,end-4))),axis image, axis off
% subplot(143),imagesc(real(h.Results.Field3D(:,:,end-2))),axis image, axis off
% subplot(144),imagesc(real(h.Results.Field3D(:,:,end-0))),axis image, axis off
cmap = cbrewer('seq','Greens',64);
cmap(cmap<0) = 0;
colormap(cmap)
% colormap viridis
%% Backward model
%% Rytov + Fourier diffraction theorem
clf('reset')
h.params_optics.verbose = false;
retAmplitude = squeeze(abs(h.Results.Field2D_trans(:, :, 1, :)));
retPhase = squeeze((angle(h.Results.Field2D_trans(:, :, 1, :))));
for j1 = 1:size(retPhase, 3)
    %%
    [~,~,~,~,a0] = PhiShiftMS(retAmplitude(:, :, j1).^2,1,0);
    retAmplitude(:, :, j1) = retAmplitude(:, :, j1) ./ sqrt(a0);
    retPhase(:, :, j1) = PhiShiftMS(Unwrap_TIE_DCT_Iter((retPhase(:, :, j1))), 1, 1);
    if h.params_optics.verbose
        ax1 = subplot(121); imagesc(retAmplitude(:, :, j1)), axis image, axis off,colorbar
        ax2 = subplot(122); imagesc(retPhase(:, :, j1)), axis image, axis off, set(ax1, 'colormap', gray), set(ax2, 'colormap', parula), set(gcf, 'color', 'w'), colorbar, drawnow
    end
end
[n_recon, ORytov,res3]=FDT_MS_simul(h,retAmplitude,retPhase,99); % FDT
n_recon = real(n_recon);
h.backward.mask = (ORytov~=0);
% n_recon(n_recon<h.params_optics.n_m)=h.params_optics.n_m;
% n_recon = NN(h,n_recon, ORytov, res3, 1);
% h.backward.outer_iteration_num = 1;
% h.backward.inner_iteration_num = 20;
% h.backward.mu = 0.01;
% error_single_scattering = [];
% for j1 = 1:500
% n_recon = TV_FISTA(h,n_recon, ORytov~=0);
% error_single_scattering = [error_single_scattering sum(abs(n_recon-h.Sample.RImap).^2,'all') ./ sum(abs(n_recon-h.params_optics.n_m).^2,'all')];
% end
% % n_recon=Nonnegative_iterations_handle(n_recon,ORytov,res3(1),res3(1),h.params_optics.lambda,h.params_optics.n_m,10, h.params_optics.use_GPU);

%% CBS backward
close all
k = h.params_optics.n_m / h.params_optics.lambda*2*pi;
h.params_optics.verbose = true;
h.params_optics.use_cuda = true;
h.backward.L1_threshold =1;
h.backward.t0 = 0.01;
h.backward.epoch_max = 500;
h.backward.t0_decay_rate = 1;
h.backward.num_pickup_scan = size(h.source.y_in,4);
h.backward.GradientMethod = 'Vanilla';
h.backward.nmax = 1.5;
h.backward.nmax_imag = 0;
% h.backward.regularizer = {'Nonnegativity',  'RImax'}; % {'Nonnegativity',  'RImax','TV'}
h.backward.regularizer = {};
h.backward.outer_iteration_num = 1;
h.backward.inner_iteration_num = 20;
h.backward.mu = 0.0025;

Results = Backward(h);



%% Side functions
function n_recon=Nonnegative_iterations_handle(n,ORytov,res3l,res3a,lambda_img,nm_img,iterations, GPU_on)
% Nonnegative_iterations_handle
% Parameter definition
if GPU_on
    n = gpuArray(single(n)); ORytov = gpuArray(single(ORytov));
end

normFact=1./(res3l^2*res3a);
ORytov_index=((abs(ORytov)==0));
%% 1. GP
n=-(2*pi*nm_img/lambda_img)^2.*(n.^2/nm_img^2-1); % RI -> scattering potential
    for jj=1:iterations
        clc, disp([num2str(jj) ' / ' num2str(iterations)])
	% GP in image space
        n((real(n)>0))=0;
    % GP in k space
        ORytov_new=fftshift(fftn(ifftshift((n))))/normFact;
        ORytov_new=ORytov_new.*ORytov_index+ORytov;
    % scattering potential -> RI
        n=fftshift(ifftn(ifftshift(ORytov_new)))*normFact;
    end

    n=nm_img*sqrt(1-n.*(lambda_img/(nm_img*2*pi))^2); % scattering potential -> RI
    n=real(n);
    n(n<nm_img)=nm_img; 
    n_recon = gather(real(n));

end
function res_img = unwrap_phase(img)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast unwrapping 2D phase image using the algorithm given in:                 %
%     M. A. Herráez, D. R. Burton, M. J. Lalor, and M. A. Gdeisat,             %
%     "Fast two-dimensional phase-unwrapping algorithm based on sorting by     %
%     reliability following a noncontinuous path", Applied Optics, Vol. 41,    %
%     Issue 35, pp. 7437-7444 (2002).                                          %
%                                                                              %
% If using this code for publication, please kindly cite the following:        %
% * M. A. Herraez, D. R. Burton, M. J. Lalor, and M. A. Gdeisat, "Fast         %
%   two-dimensional phase-unwrapping algorithm based on sorting by reliability %
%   following a noncontinuous path", Applied Optics, Vol. 41, Issue 35,        %
%   pp. 7437-7444 (2002).                                                      %
% * M. F. Kasim, "Fast 2D phase unwrapping implementation in MATLAB",          %
%   https://github.com/mfkasim91/unwrap_phase/ (2017).                         %
%                                                                              %
% Input:                                                                       %
% * img: The wrapped phase image either from -pi to pi or from 0 to 2*pi.      %
%        If there are unwanted regions, it should be filled with NaNs.         %
%                                                                              %
% Output:                                                                      %
% * res_img: The unwrapped phase with arbitrary offset.                        %
%                                                                              %
% Author:                                                                      %
%     Muhammad F. Kasim, University of Oxford (2017)                           %
%     Email: firman.kasim@gmail.com                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Ny, Nx] = size(img);

    % get the reliability
    reliability = get_reliability(img); % (Ny,Nx)

    % get the edges
    [h_edges, v_edges] = get_edges(reliability); % (Ny,Nx) and (Ny,Nx)

    % combine all edges and sort it
    edges = [h_edges(:); v_edges(:)];
    edge_bound_idx = Ny * Nx; % if i <= edge_bound_idx, it is h_edges
    [~, edge_sort_idx] = sort(edges, 'descend');

    % get the indices of pixels adjacent to the edges
    idxs1 = mod(edge_sort_idx - 1, edge_bound_idx) + 1;
    idxs2 = idxs1 + 1 + (Ny - 1) .* (edge_sort_idx <= edge_bound_idx);

    % label the group
    group = reshape([1:numel(img)], Ny*Nx, 1);
    is_grouped = zeros(Ny*Nx,1);
    group_members = cell(Ny*Nx,1);
    for i = 1:size(is_grouped,1)
        group_members{i} = i;
    end
    num_members_group = ones(Ny*Nx,1);

    % propagate the unwrapping
    res_img = img;
    num_nan = sum(isnan(edges)); % count how many nan-s and skip them
    for i = num_nan+1 : length(edge_sort_idx)
        % get the indices of the adjacent pixels
        idx1 = idxs1(i);
        idx2 = idxs2(i);

        % skip if they belong to the same group
        if (group(idx1) == group(idx2)) continue; end

        % idx1 should be ungrouped (swap if idx2 ungrouped and idx1 grouped)
        % otherwise, activate the flag all_grouped.
        % The group in idx1 must be smaller than in idx2. If initially
        % group(idx1) is larger than group(idx2), then swap it.
        all_grouped = 0;
        if is_grouped(idx1)
            if ~is_grouped(idx2)
                idxt = idx1;
                idx1 = idx2;
                idx2 = idxt;
            elseif num_members_group(group(idx1)) > num_members_group(group(idx2))
                idxt = idx1;
                idx1 = idx2;
                idx2 = idxt;
                all_grouped = 1;
            else
                all_grouped = 1;
            end
        end

        % calculate how much we should add to the idx1 and group
        dval = floor((res_img(idx2) - res_img(idx1) + pi) / (2*pi)) * 2*pi;

        % which pixel should be changed
        g1 = group(idx1);
        g2 = group(idx2);
        if all_grouped
            pix_idxs = group_members{g1};
        else
            pix_idxs = idx1;
        end

        % add the pixel value
        if dval ~= 0
            res_img(pix_idxs) = res_img(pix_idxs) + dval;
        end

        % change the group
        len_g1 = num_members_group(g1);
        len_g2 = num_members_group(g2);
        group_members{g2}(len_g2+1:len_g2+len_g1) = pix_idxs;
        group(pix_idxs) = g2; % assign the pixels to the new group
        num_members_group(g2) = num_members_group(g2) + len_g1;

        % mark idx1 and idx2 as already being grouped
        is_grouped(idx1) = 1;
        is_grouped(idx2) = 1;
    end
end
function rel = get_reliability(img)
    rel = zeros(size(img));

    % get the shifted images (N-2, N-2)
    img_im1_jm1 = img(1:end-2, 1:end-2);
    img_i_jm1   = img(2:end-1, 1:end-2);
    img_ip1_jm1 = img(3:end  , 1:end-2);
    img_im1_j   = img(1:end-2, 2:end-1);
    img_i_j     = img(2:end-1, 2:end-1);
    img_ip1_j   = img(3:end  , 2:end-1);
    img_im1_jp1 = img(1:end-2, 3:end  );
    img_i_jp1   = img(2:end-1, 3:end  );
    img_ip1_jp1 = img(3:end  , 3:end  );

    % calculate the difference
    gamma = @(x) sign(x) .* mod(abs(x), pi);
    H  = gamma(img_im1_j   - img_i_j) - gamma(img_i_j - img_ip1_j  );
    V  = gamma(img_i_jm1   - img_i_j) - gamma(img_i_j - img_i_jp1  );
    D1 = gamma(img_im1_jm1 - img_i_j) - gamma(img_i_j - img_ip1_jp1);
    D2 = gamma(img_im1_jp1 - img_i_j) - gamma(img_i_j - img_ip1_jm1);

    % calculate the second derivative
    D = sqrt(H.*H + V.*V + D1.*D1 + D2.*D2);

    % assign the reliability as 1 / D
    rel(2:end-1, 2:end-1) = 1./D;

    % assign all nan's in rel with non-nan in img to 0
    % also assign the nan's in img to nan
    rel(isnan(rel) & ~isnan(img)) = 0;
    rel(isnan(img)) = nan;
end
function [h_edges, v_edges] = get_edges(rel)
    [Ny, Nx] = size(rel);
    h_edges = [rel(1:end, 2:end) + rel(1:end, 1:end-1), nan(Ny, 1)];
    v_edges = [rel(2:end, 1:end) + rel(1:end-1, 1:end); nan(1, Nx)];
end
function [goodp2,mask, mdx, mdy,height]=PhiShiftMS(varargin)
        p2=varargin{1};
        if length(varargin)==1
            maskstyle=1;
        elseif length(varargin)==2
            maskstyle=varargin{2};
            n=1;
        elseif length(varargin)==3
            maskstyle=varargin{2};
            n=varargin{3};
        else
            maskstyle=varargin{2};
            n=varargin{3};
            input_mask=varargin{4};
        end
        
        [imY, imX]=size(p2);
%         Uimg=unwrap2(Uimg);clc;
        bsize=8;
              
            
        switch maskstyle
            case 0 %Conserve image   
                goodp2=p2;mdx=0;mdy=0;mask=0;
            case 1 %No BGmask
                mask=ones(imY,imX);
                mask(bsize:imY-bsize+1,bsize:imX-bsize+1)=0;
                switch length(varargin)
                    case 2
                        [goodp2,coefficients]=D2_LSAms(p2,1,mask);
                         mdx=coefficients(1);mdy=coefficients(2);height=coefficients(3);
                    case 3
                        [goodp2,coefficients]=D2_LSAms(p2,n,mask);
                        if n~=0
                        mdx=coefficients(n);mdy=coefficients(2*n);height=coefficients(2*n+1);
                        else
                           mdx=0;mdy=0;height=coefficients; 
                        end
                    case 4
                        mask=input_mask;
                        [goodp2,coefficients]=D2_LSAms(p2,n,mask);
                        if n~=0
                        mdx=coefficients(n);mdy=coefficients(2*n);height=coefficients(2*n+1);
                        else
                           mdx=0;mdy=0;height=coefficients; 
                        end
                end

            case 2 %selected BGmask until modified image is satisfactory
                while true
                    figure;imagesc(p2);colorbar;
                    vertices=ginput();
                    close;
                    px_s=vertices(:,1)'; py_s=vertices(:,2)';
                    mask=poly2mask(px_s,py_s,imY,imX);
                    switch nargin
                        case 3
                            [goodp2,coefficients]=D2_LSAms(p2,1,mask);
                            mdx=coefficients(1);mdy=coefficients(2);height=coefficients(3);
                        case 4
                            [goodp2,coefficients]=D2_LSAms(p2,n,mask);
                            mdx=coefficients(n);mdy=coefficients(2*n);height=coefficients(2*n+1);

                    end

                      figure;imagesc(goodp2);colorbar;
                      answer=input('satisfied?', 's');
                  
                    if strcmp(answer,'y')
                        close;
                        break;
                    end
                    close;
                end
            case 3 %selected BGmask2 until modified image is satisfactory
                while true
                    figure;imagesc(p2);colorbar;
                    vertices=ginput();
                    close;
                    px_s=vertices(:,1)'; py_s=vertices(:,2)';
                    mask=~poly2mask(px_s,py_s,imY,imX);
                    switch nargin
                        case 3
                            [goodp2,coefficients]=D2_LSAms(p2,1,mask);
                            mdx=coefficients(1);mdy=coefficients(2);height=coefficients(3);
                        case 4
                            [goodp2,coefficients]=D2_LSAms(p2,n,mask);
                            mdx=coefficients(n);mdy=coefficients(2*n);height=coefficients(2*n+1);

                    end

                      figure;imagesc(goodp2);colorbar;
                      answer=input('satisfied?', 's');
                  
                    if strcmp(answer,'y')
                        close;
                        break;
                    end
                    close;
                end
            case -1
                goodp2=p2.*(-1);mdx=0;mdy=0;mask=0;
        end
end           
function [goodp2,coefficients]=D2_LSAms(varargin)
    % Nth order ramp elimination based upon least square approximation
    % coefficient: [xn, xn-1,...x1,yn,yn-1,...y1,a0]
    % Except 0-padded map
        p2=varargin{1};
        if length(varargin)==1
            n=1;
        end
        
        switch length(varargin)
            case 2        
                n=varargin{2};
                if n~=0
                    X=zeros(imY*imX,n);Y=X;
                    for ii=1:n
                        XXX=XX.^ii;X(1:imY*imX,ii)=XXX(:);
                        YYY=YY.^ii;Y(1:imY*imX,ii)=YYY(:);
                    end
                    E=ones(imY*imX,1);AA=[X,Y,E];
                    coefficients=(AA'*AA)\(AA'*p2(:));
                else
                    coefficients=mean2(p2);
                end

            case 3
                n=varargin{2};
                mask=varargin{3};
                
                [imY, imX]=size(p2);
                [XX, YY]=meshgrid(1:imX,1:imY);

                p2mask=p2.*mask;
                p2mask=omit_outliers(p2mask);
                p2mask=p2mask(:);     
                if n~=0
                X=zeros(sum(p2mask~=0),n);Y=X;
                    for ii=1:n
                    XXX=XX.^ii;XXX=XXX(:);XXX(p2mask==0)=[];X(:,ii)=XXX;
                    YYY=YY.^ii;YYY=YYY(:);YYY(p2mask==0)=[];Y(:,ii)=YYY;
                    end
                    p2mask(p2mask==0)=[];
                    E=ones(sum(p2mask~=0),1);AA=[X,Y,E];
                    coefficients=(AA'*AA)\(AA'*p2mask);
                else
                   coefficients=sum(sum(p2mask))./sum(sum(mask)); 
                end
        end
        goodp2=p2-coefficients(end).*ones(imY,imX);
        for ii=1:n
           goodp2=goodp2-coefficients(ii).*XX.^ii;
           goodp2=goodp2-coefficients(n+ii).*YY.^ii;
        end
        
        
end
function p2mask=omit_outliers(p2mask)
    list=p2mask(:);list(list==0)=[];
    p25=prctile(list,25);p75=prctile(list,75);
    cmin=p25-1.5*(p75-p25);cmax=p75+1.5*(p75-p25);
    p2mask(p2mask<cmin)=0;
    p2mask(p2mask>cmax)=0;
end

function out_mat=regulariser_MS(reg_name,in_mat,mask, non_neg, non_pos, h)
% out_mat=regulariser_MS(reg_name,in_mat,fft_weight,lambda,non_neg,outer_itt,inner_itt,use_gpu,use_cuda)
mask = fftshift(mask);
lambda = h.backward.mu;
outer_itt = h.backward.outer_iteration_num;
inner_itt = h.backward.inner_iteration_num;
use_gpu = h.params_optics.use_GPU;
use_cuda = h.params_optics.use_cuda;
% code by Herve Hugonnet
% from the paper : 'Fast Gradient-Based Algorithms for Constrained Total Variation Image Denoising and Deblurring Problems'
% from the paper : 'Hessian Schatten-Norm Regularization for Linear Inverse Problems'
% inspired by the matlab code from Amir Beck and Marc Teboulle (c.f. https://sites.google.com/site/amirbeck314/software)

% isotropic TV is used and both hessian schatten-norm and TV ar L1
% all operation are circular symetric
% monotone version is used

% the solved problem is :
%           min_x(|fft(x).*fft_weight-fft_weight.*fft(in_mat)|^2 + 2*lambda*TV(x))
% non negativity can also be used :
%           min_x(|fft(x).*fft_weight-fft_weight.*fft(in_mat)|^2 + 2*lambda*TV(x)) + inf*(x>0)

% carefull, fft_weight is not the psf ! so the input must be already deconvolved

% reg_name -> the used regularisation : 'tv' or 'hessian'
% in_matt -> the deconvoluted input to be regularised
% fft_weight -> the regularisation strength
% lambda -> the regularisation parameter
% non_neg -> boolean of wether to use non negativity or not
% outer_itt -> the number of outer iteration
% inner_itt -> the number of inner iteration

real_input=isreal(in_mat);
max_mat = max(abs(in_mat(:)))/0.299;
in_mat = in_mat ./ max_mat;
if ~real_input
    warning('The input is complex so optimisation might be slower ! ');
end
if ~strcmp(reg_name,'tv') && ~strcmp(reg_name,'hessian')
    error('The regularisation must be either ''tv'' or ''hessian'' ')
end
if ~isempty(find(mask>1,1)) || ~isempty(find(mask<0,1))
    error('Please set the weight to value between 0 and 1')
end
if ~isempty(find(abs(in_mat)>0.3,1))
    error('For the algorithm to properly work the input must be scaled between 0.3 and -0.3');
end
if non_neg && isempty(find(in_mat<0,1))
    warning('Non negativity is used but the input has no negative value. Is the scaling correct ?');
end
if inner_itt<=0 || outer_itt<=0
    error('The number of iterations must be of at least 1 ');
end
if length(size(in_mat))>3 || length(size(in_mat))<2
    error('Only 2D and 3D matrix are supported');
end


if ~use_gpu
    use_cuda=false;
    
    in_mat=single(in_mat);
    mask=single(mask);
else
    in_mat=gpuArray(single(in_mat));
    mask=gpuArray(single(mask));
end

sz1=size(in_mat,1);
sz2=size(in_mat,2);
sz3=1;
if length(size(in_mat))==3
    sz3=size(in_mat,3);
end

mask=ifftshift(mask);

A=@(X) 1/sqrt(sz1*sz2*sz3)*fftn(X).*mask;
if real_input
    A_trans=@(X) real(sqrt(sz1*sz2*sz3)*ifftn(X.*conj(mask)));
else
    A_trans=@(X) sqrt(sz1*sz2*sz3)*ifftn(X.*conj(mask));
end

alpha=max(max(abs(mask(:)).^2));
y=A(in_mat);%the base data

if strcmp(reg_name,'tv')
    cost=@(X) sum(abs(A(X)-y).^2,'all') + 2*lambda*TV_val(X);
elseif strcmp(reg_name,'hessian')
    error('To do');
else
    error('Unrecognised regularisation name');
end

s_n=0;

t_n=0;
t_np=1;

u_n=in_mat;
x_n=in_mat;
clear in_mat;

c_n=0;
c_np=Inf;


cost_history=[];
if h.params_optics.verbose
    figure(20);
end

%start the iterations
for mm=1:outer_itt

    t_n=t_np;
    c_n=c_np;
    
    %denoising for different regularisation
    if strcmp(reg_name,'tv')
        s_n=TV_regulariser_inner(u_n-(1/alpha)*A_trans(A(u_n)-y),lambda/alpha,non_neg,non_pos,inner_itt,h,max_mat);
    elseif strcmp(reg_name,'hessian')
        error('To do');
    else
        error('Unrecognised regularisation name');
    end
    
    t_np=(1+sqrt(1+4*t_n^2))/2;
    c_np=cost(s_n);
    if c_np>c_n
        c_np=c_n;
        u_n=x_n+(t_n/t_np)*(s_n-x_n);
    else
        u_n=s_n+(t_n-1)/t_np*(s_n-x_n);
        x_n=s_n;
    end
    
    cost_history(end+1)=gather((c_np(:)));
    if h.params_optics.verbose
        plot(cost_history);title('TV cost history');drawnow;
    end
end

out_mat=x_n;

if use_gpu
    out_mat=gather(out_mat);
end

end
function in_mat=project_non_neg(in_mat,non_neg,non_pos,h,max_mat)
n_m = h.params_optics.n_m;
lambda = h.params_optics.lambda;
% compute K_3
k0=1/lambda; %[um-1, spatial frequency @ vacuum]
k =2*pi*n_m*k0; % [um-1, spatial wavenumber @ medium ]


if non_neg
    Vmin_real = k^2*(h.backward.nmin^2 / n_m^2 - 1);
    Vmin_real = Vmin_real ./ max_mat;
    in_mat(in_mat<Vmin_real)=Vmin_real;
end
if non_pos
    Vmax_real = k^2*(h.backward.nmax^2 / n_m^2 - 1);
    Vmax_real = Vmax_real ./ max_mat;
    in_mat(real(in_mat) > Vmax_real) = Vmax_real; % Set max RI
    Vmax_imag = h.backward.nmax_imag^2/n_m^2-1;
    Vmax_imag = Vmax_imag ./ max_mat;
    in_mat(imag(in_mat) > Vmax_imag) = real(in_mat(imag(in_mat) > Vmax_imag)) + 1i .* Vmax_imag; % Set max imag
end

end
function P=project_TV(P)
if length(size(P{1}))==3
    A=sqrt(max(P{1}.^2+P{2}.^2+P{3}.^2,1));
    P{1}=P{1}./A;
    P{2}=P{2}./A;
    P{3}=P{3}./A;
else
    A=sqrt(max(P{1}.^2+P{2}.^2,1));
    P{1}=P{1}./A;
    P{2}=P{2}./A;
end
end
function val=TV_val(in_mat)
P=TV_L_trans(in_mat);
if length(size(in_mat))==3
    val=sum(sqrt(abs(P{1}).^2+abs(P{2}).^2+abs(P{3}).^2),'all');%iso tv
else
    val=sum(sqrt(abs(P{1}).^2+abs(P{2}).^2),'all');%iso tv
end
end
function out_mat=TV_L(P)
out_mat=P{1}+P{2};
if length(size(out_mat))==3
    out_mat=out_mat+P{3};
end

if length(size(out_mat))==3
    out_mat=out_mat-circshift(P{1},[1 0 0]);
    out_mat=out_mat-circshift(P{2},[0 1 0]);
    out_mat=out_mat-circshift(P{3},[0 0 1]);
else
    out_mat=out_mat-circshift(P{1},[1 0]);
    out_mat=out_mat-circshift(P{2},[0 1]);
end
end
function out_mat=TV_L_trans(in_mat)
if length(size(in_mat))==3
    out_mat{1}=in_mat-circshift(in_mat,[-1 0 0]);
    out_mat{2}=in_mat-circshift(in_mat,[0 -1 0]);
    out_mat{3}=in_mat-circshift(in_mat,[0 0 -1]);
else
    out_mat{1}=in_mat-circshift(in_mat,[-1 0]);
    out_mat{2}=in_mat-circshift(in_mat,[0 -1]);
end
end
function out_mat=TV_regulariser_inner(in_mat,lambda,non_neg,non_pos,inner_itt,h,max_mat)

dim_num=length(size(in_mat));

if dim_num==3
    dividend=12; 
else
    dividend=8; 
end

P_n=0;
P_np{1}=0.*in_mat;
P_np{2}=0.*in_mat;
if dim_num==3
    P_np{3}=0.*in_mat;
end
R=P_np;

t_n=1;
t_np=1;

%start the iterations
for mm=1:inner_itt
    P_n=P_np;
    t_n=t_np;
    
    P_np=TV_L_trans(project_non_neg(in_mat-lambda*TV_L(R),non_neg,non_pos,h,max_mat));
    for kk=1:dim_num
        P_np{kk}=R{kk}+1/(dividend*lambda)*P_np{kk};
    end
    P_np=project_TV(P_np);
    
    t_np=(1+sqrt(1+4*t_n^2))/2;
    
    for kk=1:dim_num
       R{kk} = P_np{kk} + ((t_n-1)/t_np)*(P_np{kk}-P_n{kk});
    end
    
end

out_mat = project_non_neg(in_mat-lambda*TV_L(P_np),non_neg,non_pos,h,max_mat);
out_mat = out_mat .* max_mat;
end
function [phase_unwrap,N]=Unwrap_TIE_DCT_Iter(phase_wrap)   
   phi1 = unwrap_TIE(phase_wrap);
   phi1=phi1+mean2(phase_wrap)-mean2(phi1); %adjust piston
    K1=round((phi1-phase_wrap)/2/pi);  %calculate integer K
    phase_unwrap=phase_wrap+2*K1*pi; 
    residue=wrapToPi(phase_unwrap-phi1);
    phi1=phi1+unwrap_TIE(residue);
    phi1=phi1+mean2(phase_wrap)-mean2(phi1); %adjust piston
    K2=round((phi1-phase_wrap)/2/pi);  %calculate integer K
    phase_unwrap=phase_wrap+2*K2*pi; 
    residue=wrapToPi(phase_unwrap-phi1);
    N=0;
   while sum(sum(abs(K2-K1)))>0 
       K1=K2;
       phic=unwrap_TIE(residue);
     phi1=phi1+phic;
     phi1=phi1+mean2(phase_wrap)-mean2(phi1); %adjust piston
    K2=round((phi1-phase_wrap)/2/pi);  %calculate integer K
    phase_unwrap=phase_wrap+2*K2*pi; 
    residue=wrapToPi(phase_unwrap-phi1);
    N=N+1;
   end
end
function [phase_unwrap]=unwrap_TIE(phase_wrap)
      psi=exp(1i*phase_wrap);
      edx = [zeros([size(psi,1),1]), wrapToPi(diff(psi, 1, 2)), zeros([size(psi,1),1])];
      edy = [zeros([1,size(psi,2)]); wrapToPi(diff(psi, 1, 1)); zeros([1,size(psi,2)])];
       lap = diff(edx, 1, 2) + diff(edy, 1, 1); %calculate Laplacian using the finite difference
        rho=imag(conj(psi).*lap);   % calculate right hand side of Eq.(4) in the manuscript
   phase_unwrap = solvePoisson(rho); 
end
function phi = solvePoisson(rho)
    % solve the poisson equation using DCT
    dctRho = dct2(rho);
    [N, M] = size(rho);
    [I, J] = meshgrid([0:M-1], [0:N-1]);
    dctPhi = dctRho ./ 2 ./ (cos(pi*I/M) + cos(pi*J/N) - 2);
    dctPhi(1,1) = 0; % handling the inf/nan value
    % now invert to get the result
    phi = idct2(dctPhi);
end
