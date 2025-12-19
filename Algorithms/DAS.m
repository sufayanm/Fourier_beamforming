function [final_image, params] = DAS(channel_data, params)

%% Defining the initial parameters

fs    = params.sampling_frequency;
c0     = params.sound_speed ;

[N_sample, N_elements, N_waves] = size(channel_data);

if isfield(params, "modulation_frequency")
    w0 = 2*pi*params.modulation_frequency ;
else
    w0 = 0 ;
end

if isfield(params, "initial_time")
    t0 = params.initial_time ;  
else
    t0 = 0 ;
end

if isfield(params, "angle_apodization")
    angle_apodization = params.angle_apodization;
else
    angle_apodization = [];
end

[X, Z] = meshgrid(params.scan.x_axis, params.scan.z_axis) ;
X = X(:) ;
Z = Z(:) ;

N_pixels = length(X) ;

% calculate element delay
element_delay = single(sqrt((X - params.probe.x.').^2 + Z.^2)/c0);

if (abs(w0)<eps)
    channel_data = hilbert(channel_data);
end

time_axis = (0:N_sample-1)/fs + t0 ;

taper_mask = angular_apodization(X, Z, params.probe.x, angle_apodization, 'tukey', 0.25) ;


% Run Delay & Sum
tic

bf_data = das_function(channel_data, element_delay, taper_mask, N_pixels, N_waves, N_elements, w0, time_axis) ;
% bf_data = das_function_mex(channel_data, element_delay, taper_mask, N_pixels, N_waves, N_elements, w0, time_axis) ;

elapsed_time = toc();
fprintf('Elasped time for %s is %.2f seconds. \n', 'DAS', elapsed_time)

% assign phase according to 2 times the receive propagation distance
if (abs(w0) > eps)
    bf_data = bf_data.*exp(-1j*2*w0*Z/c0);
end

final_image = reshape(bf_data, [length(params.scan.z_axis), length(params.scan.x_axis)]) ;

end

