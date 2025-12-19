function [final_image, params] = DCWA(channel_data, params) 


%% Defining the initial parameters

fs    = params.sampling_frequency;
pitch = params.probe.pitch ;
c0     = params.sound_speed ;

[N_sample, N_elements, N_waves] = size(channel_data);

if isfield(params, "modulation_frequency")
    w0 = 2*pi*params.modulation_frequency ;
else
    w0 = 0 ;
end

if isfield(params, "initial_time")
    t0 = params.initial_time ;  %-4.048147907553812e-06
else
    t0 = 0 ;
end

if isfield(params, "temp_origin")
    temp_origin = params.temp_origin;
else
    temp_origin = 0;
end

if isfield(params, "angle_apodization")
    angle_apodization = params.angle_apodization;
else
    angle_apodization = [];
end

%% Setting padding values

if isfield(params, "temporal_padding")
    ntFFT = params.temporal_padding*N_sample + round(max(t0*fs)) ;
else
    params.temporal_padding = 2 ;
    ntFFT = 2*N_sample + round(t0*fs) ;
end

if rem(ntFFT,2)==1 % ntFFT must be even
    ntFFT = ntFFT+1;
    warning('Zero padding should be an even number. Temporal padding is updated to %2.0f', ntFFT)
end
fprintf('Number of samples in data - %4.0f \n', N_sample)
fprintf('Number of samples in zero-padded data - %4.0f \n', ntFFT)


if isfield(params, "spatial_padding")
    %- x-direction (same padding is used for ku and kv)
    nxFFT = ceil(params.spatial_padding*N_elements); % in order to avoid lateral edge effects
    if rem(nxFFT,2)==1 % nxFFT must be even
        nxFFT = nxFFT+1;
        warning('Zero padding should be an even number. Spatial padding is updated to %2.0f', ntFFT)
    end
else
    params.spatial_padding = 2;
    nxFFT = 2*N_elements ;
end

if rem(nxFFT-N_elements, 2)==1
    error('Please make (nxFFT-N_elements) an even number (By making N_elements even)')
end

nyFFT = nxFFT ;


if isfield(params, "z_lim")
    z_limit_factor = params.z_lim;
else
    z_limit_factor = 1.5 ;
end

%% Defining frequency co-ordinates

grd_size_x = 2 ; % this is used to define the kx and kz grid. size(kx or kz) = grd_size*size(kv or f)
grd_size_z = 2 ;

z_limit = round(z_limit_factor*params.scan.z_axis(end)/c0*fs) ;

if rem(z_limit, 2)==1
    z_limit = z_limit +1;
end

kx = single(2*pi*([0:grd_size_x/2*nxFFT-1 -grd_size_x*nxFFT/2:-1])/pitch/nxFFT) ;
old_kz = single(2*pi*([0:grd_size_z/2*z_limit-1 -grd_size_z/2*z_limit:-1])'*fs/z_limit/c0 + w0/c0) ;
kz = old_kz ;
kz(kz<0) = [] ;
[KXX, ~] = meshgrid(kx, kz) ;

kv = 2*pi*((0:nxFFT-1) - floor(nxFFT/2))/pitch/nxFFT;       % Receiving element
f_shifted = (((0:ntFFT-1) - floor(ntFFT/2)).')*2*pi*fs/ntFFT + w0 ;
[kv_mat, k_mat] = meshgrid(kv, f_shifted/c0) ;

ku = single(2*pi*reshape(((0:nyFFT-1) - floor(nyFFT/2)), 1, 1, [])/pitch/nyFFT);  % transmit element (Shifted according to MATLAB fft algorithm)
kmig = find_kmig(kz, kx, ku) ;


ku_z_mig =  sqrt(kmig.^2 - ku.^2) ;
kv_z_mig =  sqrt(kmig.^2 - (KXX - ku).^2) ;

% mult_factor = kmig./(sqrt(ku_z_mig.*kv_z_mig).*(ku_z_mig + kv_z_mig) + 1e7) ;
% mult_factor(isnan(mult_factor)|isinf(mult_factor)) = 0 ;

mult_factor = kmig./(sqrt(ku_z_mig.*kv_z_mig).*(ku_z_mig + kv_z_mig)) ;
mult_factor(imag(ku_z_mig)~=0 | imag(kv_z_mig)~=0) = 0 ;

limit = 0.0007 ;   % Limit maximum and minimum values in the multiplication factor

mult_factor(isinf(mult_factor)) = limit*sign(mult_factor(isinf(mult_factor))) ;
mult_factor(isnan(mult_factor)) = 0 ;
mult_factor(mult_factor<-limit) = -limit ;
mult_factor(mult_factor>limit) = limit;

kz_border = sqrt(abs(ku.^2 - (kx - ku).^2)) ;
mult_factor(kz<kz_border) = 0 ;

%% Temporal FFT
t_shift = round(2*temp_origin/c0*fs) ;   %% Index to which origine is to be moved

tic
channel_data = fft(channel_data, ntFFT) ;     % normal fft

%% Delay the data as per initial time and shifted origin

tmp = fftshift(f_shifted).*(t0 - (t_shift/fs));
channel_data = channel_data.*exp(-1i*tmp) ;

%% Spatial FFTs

fprintf('Taking spatial fft. \n')
channel_data = cat(2, zeros(ntFFT, (nxFFT-N_elements)/2, N_elements), channel_data, zeros(ntFFT, (nxFFT-N_elements)/2, N_elements)) ;  % To reconstruct the region outside aperture
channel_data = fft(fftshift(channel_data, 2), [], 2) ;

channel_data = cat(3, zeros(ntFFT, nxFFT, (nyFFT-N_elements)/2), channel_data, zeros(ntFFT, nxFFT, (nyFFT-N_elements)/2)) ;
channel_data = fft(fftshift(channel_data, 3), [], 3) ;

channel_data = fftshift(channel_data) ;

%% Angular apodization and Evanescent waves

if ~isempty(angle_apodization)
    taper_mask = taper_window( ku, kv, f_shifted/c0, angle_apodization, "tukey", 0.25) ;
    channel_data = channel_data.*taper_mask ;
    channel_data(abs(f_shifted/c0) < abs(ku) | abs(f_shifted/c0) < abs(kv)) = 0;    % Evanscent waves
else
    channel_data(abs(f_shifted/c0) < abs(ku) | abs(f_shifted/c0) < abs(kv)) = 0;    % Evanscent waves
end

%% Migrating the data

migSIG = zeros(length(kz), length(kx), nyFFT);  % Defining image matrix

fprintf('Interpolating %4.0f frames. \n', nyFFT)
for jj = 1:nyFFT
    migSIG(:, :, jj) = -((4*pi)^2)*interp2(kv_mat, k_mat, channel_data(:, :, jj), (KXX - ku(jj)), kmig(:, :, jj), 'linear', 0) ;
end

if t_shift ~= 0
    migSIG = migSIG.*exp(-1i*c0*kmig*t_shift/fs) ;   % Move temporal origin to its correct place
end

migSIG = migSIG.*mult_factor ;

%% Summing over ku and 2D inverse Fourier transform 
f_migSIG = sum(migSIG, 3) ;

f_migSIG = [zeros(length(old_kz)-length(kz), length(kx)); f_migSIG] ;
kz = old_kz ;
final_image = fftshift(ifft2(f_migSIG), 2) ;

%% Print time
elapsed_time = toc ;
fprintf('Elasped time for %s is %.2f seconds. \n', "Fourier-DAS", elapsed_time)

%% Defining spatial co-ordinates  
x = (0:length(kx)-1)*pitch/grd_size_x;     % If reconstructing only along the aperture, use nx instead of nxFFT
x = x - (x(end)+x(1))/2 + pitch/grd_size_x/2 ;     
z = (0:length(kz)-1)*c0/(grd_size_z*fs);    

%% Scaling with z coordinate
final_image = final_image.*z(:); 

%% Select the region defined in scan

x_scan_limits = [min(params.scan.x_axis(:)) max(params.scan.x_axis(:))] ;
z_scan_limits = [min(params.scan.z_axis(:)) max(params.scan.z_axis(:))] ;

if x_scan_limits(1)<min(x) || x_scan_limits(2)>max(x)
    x_index_limit = [1 length(x)];
else
    x_index_limit = index_finder(x_scan_limits, x(:)) ;
end
if z_scan_limits(1)<min(z) || z_scan_limits(2)>max(z)
    z_index_limit = [1 length(z)];
else
    z_index_limit = index_finder(z_scan_limits, z(:)) ;
end

x = x(x_index_limit(1):x_index_limit(2)) ;
z = z(z_index_limit(1):z_index_limit(2)) ;

final_image = final_image(z_index_limit(1):z_index_limit(2) , x_index_limit(1):x_index_limit(2)) ;

params.scan.x_axis = x ;
params.scan.z_axis = z ;

end


