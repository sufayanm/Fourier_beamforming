function [final_image, params] = Wavenumber_algorithm(channel_data, params) 

% Implements the wavenumber algorithm for ultrasound imaging.
%   
% Inputs:
%   channel_data : Time-domain RF data of size [N_samples × N_elements × N_waves]
%
%   params : Structure containing acquisition, probe, and reconstruction
%            parameters:
%       - sampling_frequency   : Sampling frequency (Hz)
%       - sound_speed          : Speed of sound (m/s)
%       - probe.pitch          : Element pitch (m)
%       - scan.x_axis          : Lateral reconstruction grid (m)
%       - scan.z_axis          : Axial reconstruction grid (m)
%       - temporal_padding     : Temporal zero-padding factor (optional)
%       - spatial_padding      : Spatial zero-padding factor  (optional)
%       - angle_apodization    : Angular apodization limit (degrees) (optional)
%       - modulation_frequency : Modulation frequecy for IQ data (optional)
%       - initial_time         : Initial time offset (s) (optional)
%       - temp_origin          : To shift temporal origin for better
%                                interpolation (m) (optional)
%       - z_lim                : Limit image reconstruction in z direction
%
% Output:
%   final_image : Reconstructed complex-valued image
%
% Implemented by Sufayan Mulani (sufayanm@ifi.uio.no)
%
%   Based on the wavenumber algorithm by Hunter et al. - “The wavenumber 
%   algorithm for full-matrix imaging using an ultrasonic array,” IEEE 
%   TUFFC, vol. 55, no. 11, pp. 2450–2462, Nov. 2008, doi: 10.1109/TUFFC.952.


%% Defining the system parameters

fs    = params.sampling_frequency;        % Sampling frequency [Hz]
pitch = params.probe.pitch;               % Element pitch [m]
c0    = params.sound_speed;               % Speed of sound [m/s]

[N_sample, N_elements, N_waves] = size(channel_data);

% Modulation frequency (for IQ data)
if isfield(params, "modulation_frequency")
    w0 = 2*pi*params.modulation_frequency ;
else
    w0 = 0 ;
end

% Initial time offset
if isfield(params, "initial_time")
    t0 = params.initial_time ;  
else
    t0 = 0 ;
end

% Temporal origin shift (for better interpolation)
if isfield(params, "temp_origin")
    temp_origin = params.temp_origin;
else
    temp_origin = 0;
end

% Angular apodization limit
if isfield(params, "angle_apodization")
    angle_apodization = params.angle_apodization;
else
    angle_apodization = [];
end

%% Setting padding values

% Temporal FFT length
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

% Spatial FFT length
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

nyFFT = nxFFT ;     % Same padding for transmit and receive dimensions

% Depth limit for migration
if isfield(params, "z_lim")
    z_limit_factor = params.z_lim;
else
    z_limit_factor = 1.5 ;
end

%% Defining frequency co-ordinates

grd_size_x = 2 ;  % Grid oversampling factor in x 
grd_size_z = 2 ;  % Grid oversampling factor in x 

% Maximum depth in samples
z_limit = round(z_limit_factor*params.scan.z_axis(end)/c0*fs) ;
if rem(z_limit, 2)==1
    z_limit = z_limit +1;
end

% Lateral and axial image wavenumbers
kx = single(2*pi*([0:grd_size_x/2*nxFFT-1 -grd_size_x*nxFFT/2:-1])/pitch/nxFFT) ;
temp_kz = single(2*pi*([0:grd_size_z/2*z_limit-1 -grd_size_z/2*z_limit:-1])'*fs/z_limit/c0 + w0/c0) ;
kz = temp_kz ;
kz(kz<0) = [] ;    % keep only positive kz during calculations
[KXX, ~] = meshgrid(kx, kz) ;

% channel data receive (kv), transmit (ku) and temporal (omega) frequencies
kv = 2*pi*((0:nxFFT-1) - floor(nxFFT/2))/pitch/nxFFT;       % Receiving element
omega = (((0:ntFFT-1) - floor(ntFFT/2)).')*2*pi*fs/ntFFT + w0 ;
[kv_mat, k_mat] = meshgrid(kv, omega/c0) ;

ku = single(2*pi*reshape(((0:nyFFT-1) - floor(nyFFT/2)), 1, 1, [])/pitch/nyFFT);  % transmit element (Shifted according to MATLAB fft algorithm)

% Migrated wavenumber
kmig = find_kmig(kz, kx, ku) ;

%% Calculate multiplication factor
mult_factor = sqrt(((kmig).^2 - ku.^2).*((kmig).^2 - (KXX-ku).^2)) ;

kz_border = sqrt(abs(ku.^2 - (kx - ku).^2)) ;
mult_factor(kz<kz_border) = 0 ;    % Remove non-physical region

%% Temporal FFT

tic
channel_data = fft(channel_data, ntFFT) ;    

%% Apply initial time and origin shift in frequency domain

t_shift = round(2*temp_origin/c0*fs) ;   
tmp = fftshift(omega).*(t0 - (t_shift/fs));
channel_data = channel_data.*exp(-1i*tmp) ;

%% Spatial FFTs (receive and transmit dimensions)

fprintf('Taking spatial fft. \n')
% Zero-pad receive aperture
channel_data = cat(2, zeros(ntFFT, (nxFFT-N_elements)/2, N_elements), channel_data, zeros(ntFFT, (nxFFT-N_elements)/2, N_elements)) ;  % To reconstruct the region outside aperture
channel_data = fft(fftshift(channel_data, 2), [], 2) ;

% Zero-pad transmit aperture
channel_data = cat(3, zeros(ntFFT, nxFFT, (nyFFT-N_elements)/2), channel_data, zeros(ntFFT, nxFFT, (nyFFT-N_elements)/2)) ;
channel_data = fft(fftshift(channel_data, 3), [], 3) ;

channel_data = fftshift(channel_data) ;

%% Angular apodization and Evanescent waves suppression

if ~isempty(angle_apodization)
    apod_mask = angular_apodization_freq( ku, kv, omega/c0, angle_apodization, "tukey", 0.25) ;
    channel_data = channel_data.*apod_mask ;
    channel_data(abs(omega/c0) < abs(ku) | abs(omega/c0) < abs(kv)) = 0;    % Remove evanscent waves
else
    channel_data(abs(omega/c0) < abs(ku) | abs(omega/c0) < abs(kv)) = 0;    % Remove evanscent waves
end

%% Migration step

migSIG = zeros(length(kz), length(kx), nyFFT);  
fprintf('Interpolating %4.0f frames. \n', nyFFT)

for jj = 1:nyFFT
    migSIG(:, :, jj) = -((4*pi)^2)*interp2(kv_mat, k_mat, channel_data(:, :, jj), (KXX - ku(jj)), kmig(:, :, jj), 'linear', 0) ;
end

% Restore temporal origin
if t_shift ~= 0
    migSIG = migSIG.*exp(-1i*c0*kmig*t_shift/fs) ;   % Move temporal origin to its correct place
end

migSIG = migSIG.*mult_factor ;

%% Sum over transmit wavenumber ku and inverse FFT
f_migSIG = sum(migSIG, 3) ;

% Restore full kz axis
f_migSIG = [zeros(length(temp_kz)-length(kz), length(kx)); f_migSIG] ;
kz = temp_kz ;

final_image = fftshift(ifft2(f_migSIG), 2) ;

%% Print time
elapsed_time = toc ;
fprintf('Elasped time for %s is %.2f seconds. \n', 'WA', elapsed_time)

%% Define spatial coordinates 
x = (0:length(kx)-1)*pitch/grd_size_x;     
x = x - (x(end)+x(1))/2 + pitch/grd_size_x/2 ;     
z = (0:length(kz)-1)*c0/(grd_size_z*fs);    

%% Crop to scan region 

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


