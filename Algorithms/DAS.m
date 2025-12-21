function [final_image, params] = DAS(channel_data, params)

% Implements the Delay-and-sum algorithm for ultrasound imaging.
%   
% Inputs:
%   channel_data : Time-domain RF data of size [N_samples × N_elements × N_waves]
%
%   params : Structure containing acquisition, probe, and reconstruction
%            parameters:
%       - sampling_frequency   : Sampling frequency (Hz)
%       - sound_speed          : Speed of sound (m/s)
%       - probe.x              : Lateral positions of array elements [m]
%       - scan.x_axis          : Lateral reconstruction grid (m)
%       - scan.z_axis          : Axial reconstruction grid (m)
%       - angle_apodization    : Angular apodization limit (degrees) (optional)
%       - modulation_frequency : Modulation frequecy for IQ data (optional)
%       - initial_time         : Initial time offset (s) (optional)
%
% Output:
%   final_image : Reconstructed image
%
% Implemented by Sufayan Mulani (sufayanm@ifi.uio.no)


%% Defining the system parameters

fs    = params.sampling_frequency;        % Sampling frequency [Hz]
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

% Angular apodization limit
if isfield(params, "angle_apodization")
    angle_apodization = params.angle_apodization;
else
    angle_apodization = [];
end

%% Create scan grid
[X, Z] = meshgrid(params.scan.x_axis, params.scan.z_axis) ;
X = X(:) ;
Z = Z(:) ;

N_pixels = length(X) ;

time_axis = (0:N_sample-1)/fs + t0 ;

%% Compute element delays and apodization values for each pixel-element pair

element_delay = single(sqrt((X - params.probe.x.').^2 + Z.^2)/c0);

apod_mask = angular_apodization_time(X, Z, params.probe.x, angle_apodization, 'tukey', 0.25) ;

%% Run Delay & Sum
if (abs(w0)<eps)
    channel_data = hilbert(channel_data);
end

tic

bf_data = das_function(channel_data, element_delay, apod_mask, N_pixels, N_waves, N_elements, w0, time_axis) ;
% bf_data = das_function_mex(channel_data, element_delay, taper_mask, N_pixels, N_waves, N_elements, w0, time_axis) ;

% assign phase according to 2 times the receive propagation distance
if (abs(w0) > eps)
    bf_data = bf_data.*exp(-1j*2*w0*Z/c0);
end

%% Print time
elapsed_time = toc();
fprintf('Elasped time for %s is %.2f seconds. \n', 'DAS', elapsed_time)

%% Reshape beamformed data into 2D grid
final_image = reshape(bf_data, [length(params.scan.z_axis), length(params.scan.x_axis)]) ;

end

