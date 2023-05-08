%% Snow/Glacier Retreival Lookup Table Generater Using ART Model
% This program generates the 3 lookup tables that are required for the ACO hyperspectral 
% retreival algorithm (1) clean snow white sky albedo per grain size,
% (2) scaled band depth per BRF/grain size, and (3) anisotropy factor. 
% These are saved as individual files because (2) and (3) are very large
% 
% If setting up a new sensor or running a new want to resample the lookup
% tables this can be done here. The refractive index of ice is stored in
% IOP_2008_ASCIItable.txt and can be resampled to any wavelengths before
% running through the rest of the program. Bascailly just update the
% wvl_specim_full input and then run. Check units- they need to be in terms
% of meters which is converted here. 
%
% hyperResample and brf_KB12  functions needs to be in path.
%% BRF and Albedo LUT

clear; clc; close all; %housekeeping

%initialize the input variables for all simulated spectra 
theta_i = 0:5:85;  % incident illumination angles ranging 0 to 85 degrees
theta_v = 0:5:85;  % viewing illumination angles ranging 0 to 85 degrees
phi = 0:10:360;     % relative azumith angle ranging 0 to 360 degrees
grain_radius=30:30:10000; %effective grain size ranging 30 to 10,000 micrometers

% load wvl_specim_full %load Specim Wavelengths
load wvl_specim_full %load Specim Wavelengths
M=0; % clean snow

refice2008 = importdata('IOP_2008_ASCIItable.txt'); %load refractive index (Warren 2008)
wvl_range=[97: 265]; % 350 to 2600 nm
wvl = refice2008(wvl_range,1)*1000; %get wavelength from optical property file

%get refractive index
n= refice2008(wvl_range,2); %real part of refractive indes
ni=refice2008(wvl_range,3); %imaginary part of refractive index

%resample refractive index to specim wavelengths
n=hyperResample(n, wvl, wvl_specim );
ni=hyperResample(ni, wvl, wvl_specim );

wvl_specim = wvl_specim*1e-9; %convert wavelength from nanometer to meter

albedo = zeros(length(grain_radius),length(wvl_specim)); %preallocate
BRF_lut = zeros(length(theta_i),length(theta_v),length(phi),length(grain_radius),length(wvl_specim)); %preallocate

%loop through each configuration to simulate the spectra
for i = 1:length(theta_i)
    for j = 1:length(theta_v)
        for k = 1:length(phi)
            for p = 1:length(grain_radius)
                [BRF albedo] = brf_KB12(theta_i(i), theta_v(j), phi(k), grain_radius(p), M, n, ni, wvl_specim);
                BRF_lut(i,j,k,p,:)= BRF;
                albedo_lut(p,:) = albedo;
                clear BRF albedo
            end
        end       
    end
      i % ticker to ensure program is progressing
end

save('BRF_lut.mat','BRF_lut','-v7.3') %LARGE FILE
save('albedo_lut.mat','albedo_lut')

clear i j k M p refice2008 wvl wvl_range n ni
%% Generate the scaled band depth LUT for the grain size retreival 

wvl_specim_nm=wvl_specim*10^9; %convert back to nm

%find the wavelength range for ice absorption feature
[left_val,left_shoulder]=min(abs(wvl_specim_nm-900));
[right_val,right_shoulder]=min(abs(wvl_specim_nm-1130));
nbnds =length(left_shoulder:right_shoulder); % get number of bands considered in continnum

%initialize arrays
CR=zeros(length(theta_i),length(theta_v),length(phi),length(grain_radius),nbnds); %preallocate
sbd_lut=zeros(length(theta_i),length(theta_v),length(phi),length(grain_radius)); %preallocate

%loop through each configuration to simulate the spectra
for i = 1:length(theta_i)
    for j = 1:length(theta_v)
        for k = 1:length(phi)
            for p = 1:length(grain_radius)
                    
                % find the continuum for each spectra
                spec = real(squeeze(BRF_lut(i,j,k,p,left_shoulder:right_shoulder)));
                CR(i,j,k,p,:) =  continuum_removed(spec');
                    
                % find the left and right bands (shortb and longb)
                spec = squeeze(CR(i,j,k,p,:)); %assign spec  
                [val , midpt] = min(CR(i,j,k,p,:)); 

                sbd_lut(i,j,k,p) = 1-val; %save scaled band depth

            end
        end
    end
      i % ticker to ensure program is progressing
end

save('sbd_lut.mat','sbd_lut','-v7.3') %LARGE FILE
clear i j k left_val left_shoulder midpt nbnds p right_val right_shoulder val CR spec
clear sbd_lut % to save memory, not needed for anisotropy

%% Generate the anisotropy factor LUT
% NOTE: for very small reflectance values in the BRF or albeo sometimes
% will result in extremely large anisotropy factors. In these cases the
% anisotropy factor is simply set to 1. This will introduce very small
% error because the reflectance values and solar flux in these regions are
% below 0.01%

c = zeros(length(theta_i),length(theta_v),length(phi),length(grain_radius),length(wvl_specim)); %preallocate anisotropy factor

for i = 1:length(theta_i)
    for j = 1:length(theta_v)
        for k = 1:length(phi)
            for p = 1:length(grain_radius)
                for n =1:length(wvl_specim)
                    if albedo_lut(p,n)<0.001 || BRF_lut(i,j,k,p,n) <0.01 %remove large anistorpy factors due to very small numbers being divided
                        c(i,j,k,p,n) = 1; %dont scale and just use the measured BRF
                    else
                        c(i,j,k,p,n) = albedo_lut(p,n)' ./ squeeze(BRF_lut(i,j,k,p,n));
                    end
                end
            end
        end
    end
      i % ticker to ensure program is progressing
end

save('anisotropy_factor_lut.mat','c','-v7.3') %LARGE FILE