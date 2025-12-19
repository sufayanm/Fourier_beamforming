function bf_data = das_function(channel_data, element_delay, taper_mask, N_pixels, N_waves, N_elements, w0, time_axis)


% Allocate memory
bf_data=complex(zeros([N_pixels, 1], 'single'));

% transmit loop
for n_tx = 1:N_waves
    
    % receive loop
    for n_rx = 1:N_elements

        apodization = taper_mask(:,n_rx, n_tx) ;
        delay = element_delay(:, n_rx) + element_delay(:, n_tx) ;

        % beamformed signal
        temp = apodization.*interp1(time_axis, channel_data(:, n_rx, n_tx, :), delay, 'linear', 0) ;

        % apply phase correction factor to IQ data
        if(abs(w0)>eps)
            temp = exp(1i.*w0*delay).*temp ;
        end
        bf_data = bf_data + temp;
    end
end