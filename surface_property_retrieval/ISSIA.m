%% Airborne Coastal Observatory: Hyperspectral Retreivals for Mountain Snow and Glaciers 

% This code produces three surface property retreivals on a per flight line 
% basis. These retreivals are broadband albedo, optical grain radius, and radiative 
% forcing by light absorbing particles (LAPs)
% 
% Code written by Christopher Donahue for the Hakai/UNBC  Airborne Coastal Observatory. 
% For questions or corrections email donahue.christopher@outlook.com
% 
% The methodology is described in detail in "Bridging the gap between 
% airborne and spaceborne imaging spectroscopy for mountain glacier surface 
% property retrievals" under review in Remote Sensing of Environment.

% All uses of this code or derivative products from this code should use the above
% citation
% 
% This program requres the hyperspectral and statistics MATLAB toolboxes to 
% run. Additionally requires the continuum_removed function (packaged with
% program), which should be in the path.
% 
% To run the program call the ISSIA function and pass the folder that that
% contains the following ATCOR-4 output files for each flight line for processing

% ISSIA(selpath,zone)
%selpath is the full path of the folder that contains all of the ATCOR-4
%output files used for processing. The following files should be in the
%folder for each flight line.
% 1) The ATCOR-4 .inn file (this pulls the solar zenith angle for each flight line)
% 2) The atmospherically corrected reflectance file (atm.dat)
% 3) The modeled global solar flux at the ground (eglo.dat)
% 4) The terrain slope map from ATCOR-4 (slp.dat)
% 5) The terrain aspect map from ATCOR-4 (asp.dat)

%zone is the UTM zone of the flight lines being processed. For example,
%Place Glacier flight lines would use '10N'. 

%OUTPUT
% A separate geotiff is saved for each retreival in the same directory as the in put files. 
% Pixels with retreival value are stored with a nan.

% **NOTE**: This code assumes that the atm, eglo, slp, and asp files all have the 
% same spatial extent and resolution. This is acheived by performing a spatial 
% subset for each output. Spatial subsetting was required anyways since the file 
% sizes where too large to proceess all at once. 

% Three lookup tables are needed to run the code (included). All the lookup
% tables are setup for the AisaFENIX spectrometer run in binning mode to
% that outputs 451 bands from 380-2500 nm. If a different instrument of
% spectral resolution is used then these look up tables (except scaled band
% depth) need to be resampled.
% 1) The scaled band depth LUT contains the depth of the absoption feature for 
%    grain radius ranging from 30 to 5000 micrometers. 
%    It is a 4D array: [illumination angle, viewing angle, relative aziumuth, grain radius]
%     indexing for each of the parameters is included here
% 
% 2) The anisotropy factor LUT (large file size) - Contains the spectral
%    anisotorpy factor for various local illumination angles, viewing
%    angles, relative aziumuth, and grain size. It contains a 5D array (c)
%    c = [illumination angle, viewing angle, relative aziumuth, grain radius, wavelength]
%    Indexing for the input varibles are same as the sbd_LUT above, though not included here
%    NOTE: This is reloaded for every flight line as is a cause for the slow speed of the program running
%    Performance could be improved here, but my computer ran out of storage so it needed to be cleared
%
% 3) The albedo_LUT is a 2D array containing white sky albedo for each
%    modeled grain size
%    albedo_lut = [grains_radius, wavelength]

% Outside of this program, the flight lines can be mosaiced together and I
% have chosen to use an average (mean)  retreival value for overlapping flight lines which
% tends to smooth out the sometimes noticible flight line artifacts. 

%For reference a single Place Glacier flight line at 1 m spatial resolution
%takes about 40 minutes on personal laptop

function ISSIA(selpath,zone) %selpath is the full folder path

% clear; clc; close all; %housekeeping - clear previous inputs

% Prompt user to select phase directory that contains flight line files
% atm, slp, asp, eglo, inn files needed for each flight line
% selpath = uigetdir(path,'Select Directory')

% Get a list of all flight lines in the folder using the inn files
filePattern = fullfile(selpath, '*.inn'); % Load all .inn data files to struct
datafiles = dir(filePattern);

%store path to all files needed for processing
for i= 1:length(datafiles)
    flight_root{i} = datafiles(i).name(1:end-4);
    dataFullName(i).inn = fullfile(selpath,  [flight_root{i} '.inn']);
    dataFullName(i).atm = fullfile(selpath,  [flight_root{i} '_atm.dat']);
    dataFullName(i).eglo = fullfile(selpath, [flight_root{i} '_eglo.dat']);
    dataFullName(i).slp = fullfile(selpath,  [flight_root{i} '_slp.dat']);
    dataFullName(i).asp = fullfile(selpath,  [flight_root{i} '_asp.dat']);
end

%get solar zenith and azimuth angles from each inn file
for i= 1:length(datafiles)
    fid = fopen(dataFullName(i).inn);
    data = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
    fclose(fid);
    theta_i_obs(i) = round(str2num(data{1, 1}{19, 1}(1:4)));  %Solar zenith angle [degrees] 
    phi_i_obs(i) = round(str2num(data{1, 1}{19, 1}(5:13))); %Solar azimuth angle [degrees]-measured from 0 deg at true north clockwise
end

%% Get coordinate reference system code 
if isequal(zone,'8N')
    coordRefSysCode = 32608;
elseif isequal(zone,'9N')
    coordRefSysCode = 32609;
elseif isequal(zone,'10N')
    coordRefSysCode = 32610;
elseif isequal(zone,'11N')
    coordRefSysCode = 32611;
elseif isequal(zone,'12N')
    coordRefSysCode = 32612;
else
    error('UTM zone not found');
end

%% Start to loop through flight lines

for p=1:length(datafiles) %loop through all flight lines (MASTER LOOP)


    %Load flight line files
    atm = hypercube(dataFullName(p).atm);
    slp = hypercube(dataFullName(p).slp,1);
    asp = hypercube(dataFullName(p).asp,1);

    %Modify data inputs to work in Matlab more efficiently 
    cube=atm.DataCube; %store datacube as double instead of int16
    wvl_specim = double(atm.Wavelength); % (nanometers) store wvl as double instead of int16
    aspect = double(asp.DataCube); %store datacube as double instead of int16
    slope  = double(slp.DataCube);  %store datacube as double instead of int16
    clear atm slp asp 
    
    % Apply spectral smoothing to the hyperspectral datacube
    cube=smoothdata(cube,3,'sgolay',10); %get the wiggles out

    % Normalized Difference Snow Index (NDSI) to create snow mask
    [numrow , numcol, numband] = size(cube); %get size of hypercube 
    image.NDSI = zeros(numrow,numcol); %preallocate mask array

    % Find the wavelength index for the VIS (600nm) and NIR (1500 nm) bands
    % Wavelength selection per Painter et al 2013
    [~ ,NDSI_vis]=min(abs(wvl_specim-600));
    [~ ,NDSI_swir]=min(abs(wvl_specim-1500));
    NDSIval=0.87; %NDSI threshold for snow/ice

    %loop through datacube and assign NDSI < 0.87 to nan
    for i = 1:numrow
        for j = 1:numcol
        
            image.NDSI(i,j) = (cube(i,j,NDSI_vis) - cube(i,j,NDSI_swir)) / (cube(i,j,NDSI_vis) + cube(i,j,NDSI_swir));
        
            if image.NDSI(i,j) < NDSIval || cube(i,j,NDSI_swir) == 0
                image.NDSI(i,j) = nan; %set non-snow pixels to NaN
            end
    
        end
    end

    %shadow mask
    image.shadow = zeros(numrow,numcol); %preallocate mask array

    %The shadows cause an increase at short wavelengths - upward 'hoooking' shape
    %The ratio of band 1 and 560 finds the shadows
    shadow_b1 = 1 ;
    [~ , shadow_b2]=min(abs(wvl_specim-560));
    shadowval=1.0; %threshold for shadows

    for i = 1:numrow
        for j = 1:numcol
        
            image.shadow(i,j) = cube(i,j,shadow_b1) / cube(i,j,shadow_b2) ;
        
            if image.shadow(i,j) > shadowval 
                image.NDSI(i,j) = nan; %set non-snow pixels to NaN
            end
    
        end
    end

    clear NDSI_vis NDSI_swir i j
    
    %Calculate the continuum removed depth of the ice absorption feature centered at 1030 nm
    image.depth = zeros(numrow, numcol); %Preallocate arrays the mosaic 

    %find the wavelength index for the left and right shoulder of the absorption feature
    [~ ,left_shoulder]=min(abs(wvl_specim-830));
    [~ ,right_shoulder]=min(abs(wvl_specim-1130));
    image.CR=zeros(numrow,numcol,length(left_shoulder:right_shoulder)); %preallocate Continuum Removed spectra

    %Loop through hypercube spatially
    for i = 1:numrow
        for j = 1:numcol

            if isnan(image.NDSI(i,j))
                image.depth(i,j)=nan; %if non-snow pixel assign NaN to grain size map
   
            else % find the continuum for each spectra
                spec = squeeze(cube(i,j,left_shoulder:right_shoulder));
                
                try
                    image.CR(i,j,:) =  continuum_removed(spec'); %store continuum removed in image.CR
                catch
                    image.CR(i,j,:) = nan(length(spec),1); %if CR function fails store spectra as NAN because something wrong with spectra
                end
          
            if isnan(image.CR(i,j,1))  %if the continuum removed outputs nan then update NDSI and depth to nans
                image.NDSI(i,j) = nan;
                image.depth(i,j)=nan;
            else %otherwise calculate the depth and store it
                [val, ~] = min(squeeze(image.CR(i,j,:))); %find the midpoint which is the min value
                image.depth(i,j) = 1-val; %store scaled band depth
            end
            end
        end
    end

    clear spec a b c d midpt left_shoulder right_shoulder val 
    image.CR = []; %clear variables for memeory

    %Determine the effective/local illumination, viewing, and azimuth angles for each pixel 
    load sbd_LUT.mat %load scaled band depth LUT

    % Viewing from the airplane is nadir
    theta_v_obs = 0; %veiwing angle [degrees]
    phi_v_obs= phi_i_obs(p); %relative azimuth [degrees]

    %preallocate the local illumination/viewing angle geometry arrays
    angle_eff.theta_i = zeros(numrow,numcol);
    angle_eff.theta_v = zeros(numrow,numcol);
    angle_eff.raa = zeros(numrow,numcol);

    theta_i_idx = zeros(numrow,numcol); %preallocate index for illumination angle
    theta_v_idx = zeros(numrow,numcol); %preallocate index for viewing angle
    raa_idx = zeros(numrow,numcol);     %preallocate index for relative azimuth angle

    [demrow, demcol] = size(slope); %get size of dem file in case its slightly smaller

    %compute the local effective angles and store in struct
    for i = 1:numrow
        for j = 1:numcol
            
            if i>demrow || j > demcol %catches error for slight mismatch in DEM and reflectance arrays
                warning('DEM and reflectance arrays have slight size mismatch')
                image.NDSI(i,j) = nan;
                theta_i_idx(i,j) = nan; 
                theta_v_idx(i,j) = nan;
                raa_idx(i,j) = nan;
            else
            [theta_i_eff, theta_v_eff, raa_eff] =  local_viewing_angle(theta_i_obs(p), phi_i_obs(p), theta_v_obs, phi_v_obs, slope(i,j), aspect(i,j));
            angle_eff.theta_i(i,j) = round(theta_i_eff);
            angle_eff.theta_v(i,j) = round(theta_v_eff);
            angle_eff.raa(i,j) = round(raa_eff); %rounding to nearest ten for coarse version

                %if angles are over 85 set to nan - too much uncertainty
                if angle_eff.theta_i(i,j) > 85 || angle_eff.theta_v(i,j) > 85 %if oblique angles add to mask
                    image.NDSI(i,j) = nan;
                    theta_i_idx(i,j) = nan; 
                    theta_v_idx(i,j) = nan;
                    raa_idx(i,j) = nan;
 
                else
                    [~ , theta_i_holder] = min( abs( theta_i - angle_eff.theta_i(i,j) ) );
                    [~ , theta_v_holder] = min( abs( theta_v - angle_eff.theta_v(i,j) ) );
                    [~ , raa_holder]     = min( abs( phi - angle_eff.raa(i,j) ) );

                    theta_i_idx(i,j) = theta_i_holder; %store for later use in albedo section
                    theta_v_idx(i,j) = theta_v_holder;
                    raa_idx(i,j) = raa_holder;
                end
            end
        end
    end

    clear slope aspect theta_i_holder theta_v_holder raa_holder i j theta_i_eff theta_v_eff raa_eff
    
    %Assign grain size based on scaled band depth LUT
    image.gs = zeros(numrow,numcol);    %preallocate grain size arrary

    for i = 1:numrow
        for j = 1:numcol
        
            if isnan(angle_eff.theta_i(i,j)) || isnan(angle_eff.theta_v(i,j)) || isnan(angle_eff.raa(i,j))
                image.gs(i,j) = nan;
            elseif isnan(image.NDSI(i,j))
                image.gs(i,j)=nan;
            else
                sbd_gs = squeeze(sbd_lut(theta_i_idx(i,j),theta_v_idx(i,j), raa_idx(i,j), :));
                image.gs(i,j) = interp1(sbd_gs, grain_radius, image.depth(i,j),'pchip'); 

                if image.gs(i,j)<30 || image.gs(i,j)>5000 %outside of model range
                    image.gs(i,j) = nan;
                end
            end
        end
    end

    clear sbd_gs sbd_lut i j 
    
    %Convert measured HCRF to albedo using anisotropy factor and calculate broadband albedo
    hcube=hypercube(dataFullName(p).eglo); %load solar flux map
    sflux=double(hcube.DataCube); %store datacube as double instead of int16
    clear hcube
    wvl_specim_um=wvl_specim/1000; %convert wavelenth to micrometers

    load anisotropy_factor_LUT.mat %load anisotropy lut

    cube=cube/10000; %turn hypercube into decimal 
    image.albedo=zeros(numrow,numcol,numband); %preallocate spectral albedo
    image.bbalbedo=zeros(numrow,numcol); %preallocate broadband albedo

    for i=1:numrow
        for j=1:numcol
        
            if theta_i_idx(i,j) == 0 || theta_v_idx(i,j) == 0 || raa_idx(i,j)==0 || isnan(image.NDSI(i,j)) %same as mask
                image.albedo(i,j,:)=nan(length(wvl_specim),1);
                image.bbalbedo(i,j)=nan;
            else
                % spectral albedo
                [~ , gs_idx] = min( abs( grain_radius - image.gs(i,j) ) ); %find index for grain size
                image.albedo(i,j,:) = squeeze(cube(i,j,:)) .* squeeze(c(theta_i_idx(i,j), theta_v_idx(i,j),raa_idx(i,j),gs_idx,:));
   
                %broadband albedo
                ref = trapz(wvl_specim_um,squeeze(sflux(i,j,:)).* squeeze(image.albedo(i,j,:))); %reflected energy
                tot = trapz(wvl_specim_um,squeeze(sflux(i,j,:))); %total engergy
                image.bbalbedo(i,j) = round(ref/tot,3);
            end

        end
    end

    clear gs_idx ref tot c
    
    %Calculate radiative forcing due to light absorbing particles (LAPs)

    load albedo_LUT.mat %load clean snow albedo vs grain size LUT

    %Define wavelenth range for LAP radiative forcing (1000 nm) per Painter %2013 et al.
    [~ ,right_idx]=min(abs(wvl_specim-1000));

    image.rf = zeros(numrow,numcol); %preallocate radiative forcing

    for i=1:numrow
        for j=1:numcol
        
            if isnan(image.gs(i,j)) %if non value apply nan
                image.rf(i,j)=nan;
            else
                [gs_val , gs_idx] = min( abs( grain_radius - image.gs(i,j) ) ); %find index for grain size
                albedo_clean = albedo_lut(gs_idx,:)'; %find clean snow albedo spectra
                rf_diff = zeros(right_idx,1); %initialize array to store differences
        
                for k=1:right_idx %loop through wavelenths to 350 to 1000 nm
                    
                    rf_diff(k) = (albedo_clean(k) - image.albedo(i,j,k)) ;

                    if rf_diff(k) < 0
                        rf_diff(k) = 0; %cant have negative radiative forcing componenets.
                    end
                end
        
                image.rf(i,j) = 0.1*(trapz(wvl_specim_um(1:right_idx),rf_diff .* squeeze(sflux(i,j,1:right_idx)) )); %calculate RF
            
            end

        end
    end

    clear right_idx rf_diff albedo_lut

    %Write Geotiff

    clear cube %need to open up some memory

    %read georaster to get georeference
    [A,R] = readgeoraster(dataFullName(p).atm);

    %write grain size
    geotiffwrite([selpath '\' flight_root{p} '_HSI_grainsize'],image.gs,R, 'CoordRefSysCode',coordRefSysCode);

    %write raditive forcing
    geotiffwrite([selpath '\' flight_root{p} '_HSI_radiativeforcing'],image.rf,R, 'CoordRefSysCode',coordRefSysCode);

    %write albedo
    geotiffwrite([selpath '\' flight_root{p} '_HSI_albedo.tif'],image.bbalbedo,R, 'CoordRefSysCode',coordRefSysCode);

end

end %end function