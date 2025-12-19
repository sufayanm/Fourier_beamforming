function taper_mask = taper_window(ku, kv, k, angle_apodization, taper_type, taper_length)

if nargin<4
    taper_type= '' ;
end

if strcmp(taper_type, "linear")
    if nargin<5
        taper_length = angle_apodization/10 ;
    end

    angle_u = rad2deg(asin(ku./k)) ;
    angle_v = rad2deg(asin(kv./k)) ;

    theta_cutoff = angle_apodization;
    theta_taper = angle_apodization - taper_length;     % Start tapering here

    taper_mask_u = ones(size(angle_u));  % Start with all ones
    taper_mask_u(abs(angle_u) > theta_cutoff) = 0;

    taper_mask_v = ones(size(angle_v));  % Start with all ones
    taper_mask_v(abs(angle_v) > theta_cutoff) = 0;

    taper_zone = (abs(angle_u) > theta_taper) & (abs(angle_u) <= theta_cutoff);
    taper_mask_u(taper_zone) = (theta_cutoff - abs(angle_u(taper_zone))) / ...
        (taper_length);

    taper_zone = (abs(angle_v) > theta_taper) & (abs(angle_v) <= theta_cutoff);
    taper_mask_v(taper_zone) = (theta_cutoff - abs(angle_v(taper_zone))) / ...
        (taper_length);

elseif strcmp(taper_type, "tukey")
    if nargin<5
        taper_length = 0.25 ;
    end

    tan_u = ku./sqrt(k.^2 - ku.^2) ;
    tan_v = kv./sqrt(k.^2 - kv.^2) ;

    tan_u(isnan(tan_u)) = 0 ;
    tan_v(isnan(tan_v)) = 0 ;

    f_number = cot(deg2rad(angle_apodization))/2 ;

    ratio_u = abs(tan_u.*f_number) ;
    ratio_v = abs(tan_v.*f_number) ;

    taper_mask_u = (ratio_u<=(1/2*(1-taper_length))) + (ratio_u>(1/2*(1-taper_length))).*(ratio_u<(1/2)).*0.5.*(1+cos(2*pi/taper_length*(ratio_u-taper_length/2-1/2)));
    taper_mask_v = (ratio_v<=(1/2*(1-taper_length))) + (ratio_v>(1/2*(1-taper_length))).*(ratio_v<(1/2)).*0.5.*(1+cos(2*pi/taper_length*(ratio_v-taper_length/2-1/2)));

else
    warning("Apodization window is not defined or does not exists. \n")
end
taper_mask = taper_mask_u.*taper_mask_v ;

taper_mask(abs(k)<abs(ku) | abs(k)<abs(kv)) = 0 ;
end