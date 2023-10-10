function imRef = high_pass_filter(h,imRef, diameter_um)
% h = TFM; imRef = FL.Ref; imMoving = FL.Moving; ROI_mask = ROI_mask_FL; num_features = 10; num_candidates = 100;
%% Noise filtering (https://www.sciencedirect.com/science/article/pii/S0021979796902179)

if nargin == 2
    diameter_um = [1 1 1];
elseif nargin == 3
    if length(diameter_um) == 1
        diameter_um = ones(1,3)*diameter_um;
    end
end

if size(imRef,3) > 1
    cfilter = zeros(round(diameter_um./h.parameters.resolution),'single');
    
    ii = 1:size(cfilter,1); ii = reshape(ii - mean(ii),[],1,1);
    jj = 1:size(cfilter,2); jj = reshape(jj - mean(jj),1,[],1);
    kk= 1:size(cfilter,3); kk = reshape(kk - mean(kk),1,1,[]);
    
    cfilter = exp(-(ii.^2+jj.^2+kk.^2)/4);
    B = sum(exp(-(ii.^2)/4)).*sum(exp(-(jj.^2)/4)).*sum(exp(-(kk.^2)/4));
    K0 = 1/B*sum(exp(-(ii.^2)/2)).*sum(exp(-(jj.^2)/2)).*sum(exp(-(kk.^2)/2))  -1 / length(cfilter(:));
    
    cfilter = 1/K0 .* (cfilter/B - 1 / length(cfilter(:)));
    
    imRef = max(imfilter(imRef, cfilter),0);
else
    cfilter = zeros(round(1./h.parameters.resolution(1:2)),'single');
    
    ii = 1:size(cfilter,1); ii = reshape(ii - mean(ii),[],1,1);
    jj = 1:size(cfilter,2); jj = reshape(jj - mean(jj),1,[],1);
    
    cfilter = exp(-(ii.^2+jj.^2)/4);
    B = sum(exp(-(ii.^2)/4)).*sum(exp(-(jj.^2)/4));
    K0 = 1/B*sum(exp(-(ii.^2)/2)).*sum(exp(-(jj.^2)/2))  -1 / length(cfilter(:));
    
    cfilter = 1/K0 .* (cfilter/B - 1 / length(cfilter(:)));
    
    imRef = max(imfilter(imRef, cfilter),0);
end
end