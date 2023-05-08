 ISSIA
 Imaging Spectrometer - Snow and Ice Algorithm
 
 This code produces three surface property retreivals on a per flight line 
 basis. These retreivals are broadband albedo, optical grain radius, and radiative 
 forcing by light absorbing particles (LAPs)
 
 Code written by Christopher Donahue for the Hakai/UNBC  Airborne Coastal Observatory. 
 For questions or corrections email donahue.christopher@outlook.com
 
 The methodology is described in detail in "Bridging the gap between 
 airborne and spaceborne imaging spectroscopy for mountain glacier surface 
 property retrievals" under review in Remote Sensing of Environment.

 All uses of this code or derivative products from this code should use the above
 citation
 
 This program requres the hyperspectral and statistics MATLAB toolboxes to 
 run. Additionally requires the continuum_removed function (packaged with
 program), which should be in the path.
 
 Running the program will prompt the user to choose a directy to pull input 
 files and save output files. The program will process each flight line in the 
 folder and the follwoing ATCOR-4 output files for each flight line need to be 
 located in the folder.
 To run the program call the ISSIA function and pass the folder that that
 contains the following ATCOR-4 output files for each flight line for processing

 INPUTS:
 1) The ATCOR-4 .inn file (this pulls the solar zenith angle for each flight line)
 2) The atmospherically corrected reflectance file (atm.dat)
 3) The modeled global solar flux at the ground (eglo.dat)
 4) The terrain slope map from ATCOR-4 (slp.dat)
 5) The terrain aspect map from ATCOR-4 (asp.dat)

OUTPUT
 A separate geotiff is saved for each retreival in the same directory as the in put files. 
 Pixels with retreival value are stored with a nan.
 NOTE: The Coordinate reference system code is a required input and is
 currently set to WGS 84 / UTM Zone 10N for Place Glacier. If processing
 data elsewhere update coordRefSysCode.

 **NOTE**: This code assumes that the atm, eglo, slp, and asp files all have the 
 same spatial extent and resolution. This is acheived by performing a spatial 
 subset for each output. Spatial subsetting was required anyways since the file 
 sizes where too large to proceess all at once. 

 Three lookup tables are needed to run the code. These can be downloaded
 from the following data repositry:  10.6084/m9.figshare.22777184
 
 All the lookup tables are setup for the AisaFENIX spectrometer run in binning mode to
 that outputs 451 bands from 380-2500 nm. If a different instrument of
 spectral resolution is used then these look up tables (except scaled band
 depth) need to be resampled.
 1) The scaled band depth LUT contains the depth of the absoption feature for 
    grain radius ranging from 30 to 5000 micrometers. 
    It is a 4D array: [illumination angle, viewing angle, relative aziumuth, grain radius]
     indexing for each of the parameters is included here
 
 2) The anisotropy factor LUT (large file size) - Contains the spectral
    anisotorpy factor for various local illumination angles, viewing
    angles, relative aziumuth, and grain size. It contains a 5D array (c)
    c = [illumination angle, viewing angle, relative aziumuth, grain radius, wavelength]
    Indexing for the input varibles are same as the sbd_LUT above, though not included here
    NOTE: This is reloaded for every flight line as is a cause for the slow speed of the program running
    Performance could be improved here, but my computer ran out of storage so it needed to be cleared

 3) The albedo_LUT is a 2D array containing white sky albedo for each
    modeled grain size
    albedo_lut = [grains_radius, wavelength]

 Outside of this program, the flight lines can be mosaiced together and I
 have chosen to use an average (mean)  retreival value for overlapping flight lines which
 tends to smooth out the sometimes noticible flight line artifacts. 

For reference a single Place Glacier flight line at 1 m spatial resolution
takes about 40 minutes on personal laptop

An example ACO dataset over Place Glacier is available for download here: 10.6084/m9.figshare.22777184
