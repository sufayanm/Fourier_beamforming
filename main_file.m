clc; clear;

addpath("Algorithms\")
addpath("Functions\")

%% Load the phantom data

load Phantom_data.mat

%% Run wavenumber algorithm

[wa_image, params_wa] = Wavenumber_algorithm(channel_data, params) ;

% Normalize the image with the maximum value
wa_image = abs(wa_image) ;
wa_image = wa_image/max(wa_image(:)) ;

%% Run DCWA

[dcwa_image, params_dcwa] = DCWA(channel_data, params) ;

% Normalize the image with the maximum value
dcwa_image = abs(dcwa_image) ;
dcwa_image = dcwa_image/max(dcwa_image(:)) ;

%% Run DAS

[das_image, params_das] = DAS(channel_data, params) ;

% Normalize the image with the maximum value
das_image = abs(das_image) ;
das_image = das_image/max(das_image(:)) ;

%% Display phantom images reconstructed using WA, DAS, and DCWA

Y_lim = [20 80] ;      
X_lim = [-18 18] ;

figure;
t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');

nexttile
imagesc(params_wa.scan.x_axis*1e3, params_wa.scan.z_axis*1e3, db(wa_image), [-60, 0])
colormap gray
axis equal tight
title("WA")
ylim(Y_lim)
xlim(X_lim)


nexttile
imagesc(params_das.scan.x_axis*1e3, params_das.scan.z_axis*1e3, db(das_image), [-60, 0])
colormap gray
axis equal tight
title("DAS")
ylim(Y_lim)
xlim(X_lim)
yticklabels([])

xx_label = xlabel(gca, 'x [mm]') ;


nexttile
imagesc(params_dcwa.scan.x_axis*1e3, params_dcwa.scan.z_axis*1e3, db(dcwa_image), [-60, 0])
colormap gray
axis equal tight
title("DCWA")
ylim(Y_lim)
xlim(X_lim)
yticklabels([])

clr_bar = colorbar ;
clr_bar.TickDirection = "out" ;

fontsize(15,'points')

ylabel(t1, 'y [mm]') ;
