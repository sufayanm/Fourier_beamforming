function taper_mask = angular_apodization_time(X, Z, ele_x, angle_apodization, taper_type, taper_length)

if nargin<4
    taper_type= '' ;
end

x_dist = X - ele_x(:).' ;

if strcmp(taper_type, "linear")
    if nargin<5
        taper_length = angle_apodization/10 ;
    end

    angle_rec = rad2deg(atan(x_dist./Z)) ;

    theta_cutoff = angle_apodization;
    theta_taper = angle_apodization - taper_length;     % Start tapering here

    taper_mask_rec = ones(size(angle_rec));  % Start with all ones
    taper_mask_rec(abs(angle_rec) > theta_cutoff) = 0;

    taper_zone = (abs(angle_rec) > theta_taper) & (abs(angle_rec) <= theta_cutoff);
    taper_mask_rec(taper_zone) = (theta_cutoff - abs(angle_rec(taper_zone))) / ...
        (taper_length);


elseif strcmp(taper_type, "tukey")
    if nargin<5
        taper_length = 0.25 ;
    end

    tan_rec = x_dist./Z ;
    tan_rec(isnan(tan_rec)) = 0 ;
    f_number = cot(deg2rad(angle_apodization))/2 ;
    ratio_rec = abs(tan_rec.*f_number) ;
    taper_mask_rec = (ratio_rec<=(1/2*(1-taper_length))) + (ratio_rec>(1/2*(1-taper_length))).*(ratio_rec<(1/2)).*0.5.*(1+cos(2*pi/taper_length*(ratio_rec-taper_length/2-1/2)));

else
    warning("Apodization window is not defined or does not exists. \n")
end

taper_mask_tran = reshape(taper_mask_rec, [length(X), 1, length(ele_x)]) ;
taper_mask = taper_mask_rec.*taper_mask_tran ;

end