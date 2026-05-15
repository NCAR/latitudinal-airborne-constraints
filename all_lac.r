## script to run all necessary steps in LAC processing and plotting

## First, download necessary data files 
# 1) download ATom 10-second merge data set from https://doi.org/10.3334/ORNLDAAC/1925 to LOCALDATA
# 2) download ATom ObsPack GV+ v6.1 file from https://doi.org/10.25925/20201204 to LOCALDATA
# 3) download OCO-2 v10 MIP posterior fluxes from https://www.gml.noaa.gov/ccgg/OCO2_v10mip/download.php to LOCALDATA/V10MIP, preserving experiment subdirectories
# 4) download OCO-2 v10 MIP posterior ATom concentration ObsPack file from https://www.gml.noaa.gov/ccgg/OCO2_v10mip/download.php to LOCALDATA/V10MIP
# 5) download ocean pCO2 flux products from https://doi.org/10.5281/zenodo.14639761 to LOCALDATA/GCBOCEAN
# 6) download GridFED fossil-fuel emissions from https://doi.org/10.5281/zenodo.13909046 to LOCALDATA/GCBFOSS
# 7) download the GCB2024 inversion fluxes from https://doi.org/10.18160/4R5W-VNBV to LOCALDATA
# 8) download NOAA annual mean global N2O file from https://gml.noaa.gov/webdata/ccgg/trends/n2o/n2o_annmean_gl.txt to LOCALDATA
# 9) download NOAA monthly mean MLO CO2 file from https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/txt/co2_mlo_surface-flask_1_ccgg_month.txt to LOCALDATA
# 10) download NOAA monthly mean global CO2 file from https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_gl.csv to LOCALDATA

## Files already existing in LOCALDATA include:
# 11) BLUE and LUCE land use land use change estimates provided by Clemens Schwingshackl
#     (BLUE_LUH2-GCB2024_ELUC_gridded_net_2016-2018_with-peat.nc and LUCE_LUH2-GCB2024_ELUC_gridded_net_2016-2018_with-peat.nc)
# 12) ATom locally generated Medusa merge files (MEDUSA_MERGE_ATOM*_200918.tbl)
# 13) Times for filtering locally influenced boundary layer samples in two formats (atom_filter_times.csv)


## Then run this script either from the R command line ("source('all_lac.r')") or the linux command line ("R CMD BATCH all_lac.r")

# 1) merge selected aircraft campaign variables onto ObsPack ID
setwd('OBS')
source('aircraft_obspack_merge.r')
aomerge()

# 2) merge model output onto selected aircraft campaign variables using ObsPack ID
setwd('../V10MIP')
source('model_aircraft_merge.r')

# 3) process all V10MIP fluxes, aggregating by 10 degress latitude
source('flux_proc.r') # (43 min)

# 4) read in flux_proc.r output and plot latitudinal gradients 
source('flux_v_lat.r')

# 5) read in flux_proc.r output and plot regional and component relationships
source('flux_v_flux.r')
fvf() # start=ISOdate(2015,1,1,tz='UTC'),end=ISOdate(2020,12,31,tz='UTC') ## entire v10MIP period

# 6) calculate and plots harmonic fits to observations and model output
setwd('../')
source('noaa_stl.r') # calculates smoothed deseasonalized trend at MLO for use in detrending obs and models
stlsite('mlo')
source('all_harm.r') # calls:
# read_airborne_obs.r (reads in airborne 10-sec merge product and trims for selected variables, only need to call once)
# harm_fit.r (bins and fits harmonics), calls:
# 	filter_airborne.r (filters obs or models for local influences and stratosphere)
# harm_plot.r (plots results of harm_fit.r)

# 7) read in harm_fit.r output and plot column-mean latitudinal concentration gradients 
source('plot_atom_col.r') # has to be run before flux_conc_cor.r for Fig2b

# 8) calculate correlations between individual press x lat binned annual mean concentrations and zonal-mean fluxes
source('flux_conc_cor.r') # calls:
# xsect_plot.r (plots cross-sections)

# 9) calculate correlations between press x lat regions and zonal-mean fluxes, and plots and outputs ECFC results
source('flux_conc_ecfc.r') # calls:
# ecplot_mlr.r (plots ECFCs and reports back results)

# 10) calculate zonal mean GCB ocean and fossil fluxes and uncertainties 
source('flux_calc_ocean_foss.r')

# 11) calculate zonal mean GCB land, LULUC, and river fluxes and uncertainties 
source('flux_calc_land_luluc_river.r')
