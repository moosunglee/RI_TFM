function params= BASIC_OPTICAL_PARAMETER(init_params)
params=struct;
params.size=[161 161 81];%size of the refractive index
params.wavelength=0.532;%wavelength
params.NA=1.4;%objective lens NA
params.RI_bg=1.3355;%refractive index out of the sample
params.resolution=[0.1 0.1 0.1];%resolution of one voxel
params.vector_simulation=false;%use polarised field or scalar field
params.use_abbe_sine=false;
params = orderfields(params);
if nargin==1
    params=update_struct(params,init_params);
end

end