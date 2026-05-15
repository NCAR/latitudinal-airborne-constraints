## Program to read in and process independent DGVM, LULUC, and river flux estimates, and adjustments 
# writes out output for Table S5 and Table 1

library('ncdf4')


## GCB2024 DGVMs:

## The GCB2024 DGVM fluxes are aggregated to zonal bands dividing at 30 S and 30 N, but we need numbers divided at 20 S and 20 N.
## It is not possible to spatially aggregate the gridded GCB2024 DGVM product, which has 'per area' fluxes, because you need to use each model's native land/sea mask, which is not provided.
## As of February 2026, the gridded GCB2024 (TRENDY v13) output with native land/sea masks is still not available online. 
## Instead, Mike O'Sullivan calculated 20 S / 20 N division values for us, for June 2016 - May 2018
## For this, Mike used only 17 of 20 models, exluding three that only reported annual NBP (DLEM, LPJml, and LPJ-GUESS)
## As of May 2026, the gridded GCB2024 values are now available online but lack complete information on areas and land fractions for every model to aggregate /m^2 fluxes. 
## Using best guesses, got pretty close (NET: 1.53 +/- 0.47, TRP: 1.83 +/- 0.59, SET: 0.24 +/- 0.14), but sticking with Mike's provided numbers.

## Via first email without uncertainty calculations:
#NET 1.514
#TRP 1.788
#SET 0.220
## Via second email with uncertainties (mean and SD of 20 DGVMs)
#North = 1.51 +/- 0.458
#Tropics = 1.79 +/- 0.642
#South = 0.22 +/- 0.146

nbpglb=(1.514+1.788+0.220)*-1
nbpnet=1.514*-1
nbptrp=1.788*-1
nbpset=0.220*-1

uncnbpglb=0.81 # calculated for 2016-2018 from Global_Carbon_Budget_2024_v1.0-3.xlsx downloaded from https://globalcarbonbudget.org/archive
uncnbpnet=0.458
uncnbptrp=0.642
uncnbpset=0.146

print('DGVM')
print('(GCB2024)')
print(paste('GLB:',round(nbpglb,2),'+/-',round(uncnbpglb,2)))
print(paste('NET:',round(nbpnet,2),'+/-',round(uncnbpnet,2)))
print(paste('TRP:',round(nbptrp,2),'+/-',round(uncnbptrp,2)))
print(paste('SET:',round(nbpset,2),'+/-',round(uncnbpset,2)))

# print((uncnbpnet^2+uncnbptrp^2+uncnbpset^2)^0.5)
# [1] 0.8020249
## similar to global value, so zones look to be uncorrelated


## LCE:

## Friedlingstein et al., Nature, 2026 reported an anthropogenic increase in lateral transfer of carbon from terrestrial organic matter to inland aquatic systems and then to the atmoshere of 0.07 +/- 0.06
# since this is not represented in DGVMs it is an adjustment to the DGVM fluxes (reduction of Sland)
# this is not the same as river fluxes, which assume a constant 0.65 from land to ocean and out of the ocean
## This adjustment is based on Tian et al., 2023: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023GB007776, which shows in Figure 7k and 7l that the increased evasion from inland waters is almost entirely N of 20 N

lceglb=0.07
lcenet=0.07
lcetrp=0
lceset=0

unclceglb=0.06
unclcenet=0.06
unclcetrp=0
unclceset=0

print('LCE')
print('(from Friedlingstein et al. (Nature, 2026) and Tian et al. (GBC, 2023))')
print('(based on Tian et al. (GBC, 2023) Fig. 7k-l, assume all is NET)')
print('GLB: 0.07 +/- 0.06')
print('NET: 0.07 +/- 0.06')
print('TRP: 0.00 +/- 0.00')
print('SET: 0.00 +/- 0.00')
print(paste('GLB:',round(lceglb,2),'+/-',round(unclceglb,2)))
print(paste('NET:',round(lcenet,2),'+/-',round(unclcenet,2)))
print(paste('TRP:',round(lcetrp,2),'+/-',round(unclcetrp,2)))
print(paste('SET:',round(lceset,2),'+/-',round(unclceset,2)))


## RSS:

## Friedlingstein et al., Nature, 2026 reported a decrease in S_land of 0.5 ± 0.3 over 2014-2023 after considering changes in land cover
## O'Sullivan et al., GBC, 2026 got 0.57 ± 0.20 for 2015-2024 but provided gridded (per area), regional (10 RECCAP), and global files
## to use gridded files, would need model-specific areas, which are not provided
## also do not have the definition of region number in the files
## instead using O'Sullivan et al. (2026) Figure 2a, assigning:
# North America 15% tropical, 85% NET
# South America 67% tropical, 33% SET
# Europe 100% NET
# Russia 100% NET
# Africa 95% tropical, 5% SET
# South Asia 50% tropical, 50% NET
# East Asia 100% NET
# Southeast Asia 100% tropical
# West/Central Asia 100% NET
# Australasia 10% tropical, 90% SET

## using O'Sullivan et al. (2026) Figure S1a, assigning:
africa=.120
northamer=.115
seasia=.090
southamer=.075
southasia=.055
eastasia=.045
austral=.030
russia=.025
europe=.010
westasia=.005
rssset=0.33*southamer+0.05*africa+0.9*austral
rsstrp=0.15*northamer+0.67*southamer+0.95*africa+0.5*southasia+seasia+0.1*austral
rssnet=0.85*northamer+europe+russia+0.5*southasia+eastasia+westasia
rssglb=rssset+rsstrp+rssnet # 0.57

## for uncertainty use sqrt(flux fraction)
## since all fluxes are positive, do not need to use abs values
uncrssglb=0.3

print('RSS')
print('(from O\'Sullivan et al. (GBC 2026) Fig. S1a)')
print(paste('GLB:',round(rssglb,2),'+/-',round(uncrssglb,2)))
uncrssnet=uncrssglb*(rssnet/rssglb)^0.5
print(paste('NET:',round(rssnet,2),'+/-',round(uncrssnet,2)))
uncrsstrp=uncrssglb*(rsstrp/rssglb)^0.5
print(paste('TRP:',round(rsstrp,2),'+/-',round(uncrsstrp,2)))
uncrssset=uncrssglb*(rssset/rssglb)^0.5
print(paste('SET:',round(rssset,2),'+/-',round(uncrssset,2)))
print((uncrssnet^2+uncrsstrp^2+uncrssset^2)^0.5)
uncrsssettrp=uncrssglb*((rssset+rsstrp)/rssglb)^0.5
print(paste('SET+T:',round(rssset+rsstrp,2),'+/-',round(uncrsssettrp,2)))


## LULUC:

## Clemens Schwingshackl provided (June 8, 2025) gridded annual means for 2016, 2017, and 2018 for BLUE and LUCE (do not exist monthly and for other 2 only gridded by country)
## peat drainage and peat fires were added from other sources (see Friedlingstein et al., 2024)
latdiv=c(-20,20)

bluefile='LOCALDATA/BLUE_LUH2-GCB2024_ELUC_gridded_net_2016-2018_with-peat.nc'
lucefile='LOCALDATA/LUCE_LUH2-GCB2024_ELUC_gridded_net_2016-2018_with-peat.nc'

bluenc=nc_open(bluefile)
lat=ncvar_get(bluenc,'lat') # same for both
lon=ncvar_get(bluenc,'lon') # same for both
area=ncvar_get(bluenc,'cell_area') # only given for BLUE, assume same
lucenc=nc_open(lucefile)
#dimensions:
#        lat = 720 ;
#        lon = 1440 ;
#        time = 3 ;
#variables:
#        float ELUC(time, lat, lon) ;
#                ELUC:standard_name = "Carbon emissions from land use change" ;
#                ELUC:units = "t C/ha/yr" ;
#        float cell_area(lat, lon) ;
#                cell_area:units = "ha" ;

blue=ncvar_get(bluenc,'ELUC')
blue=apply(blue,c(1,2),mean)
blue=blue*area # tC/yr

luce=ncvar_get(lucenc,'ELUC')
luce=apply(luce,c(1,2),mean)
luce=luce*area # tC/yr

nc_close(bluenc); nc_close(lucenc)

blueglb=sum(blue,na.rm=T)/1E9
bluenet=sum(blue[,lat>latdiv[2]],na.rm=T)/1E9
bluetrp=sum(blue[,lat<latdiv[2]&lat>latdiv[1]],na.rm=T)/1E9
blueset=sum(blue[,lat<latdiv[1]],na.rm=T)/1E9

luceglb=sum(luce,na.rm=T)/1E9
lucenet=sum(luce[,lat>latdiv[2]],na.rm=T)/1E9
lucetrp=sum(luce[,lat<latdiv[2]&lat>latdiv[1]],na.rm=T)/1E9
luceset=sum(luce[,lat<latdiv[1]],na.rm=T)/1E9

lulucglb=mean(c(blueglb,luceglb))
lulucnet=mean(c(bluenet,lucenet))
luluctrp=mean(c(bluetrp,lucetrp))
lulucset=mean(c(blueset,luceset))

#> c(lulucset,luluctrp,lulucnet)
#[1] 0.07991745 0.99944756 0.17875633

## per email from Clemens:
## "the LULUC emissions from BLUE and LUCE are 0.12 PgC/yr higher in 2016-2018 than the GCB2024 average, as the other two models used in GCB2024 (OSCAR and H&C23) have lower emissions in that period (1.01 PgC/yr)."
#> lulucglb
#[1] 1.258121
# mean would be (1.258121+1.01)/2 = 1.134, actual mean from spreadsheet is 1.136, and 1.258-1.136 = 0.122

## adjust flux proportionally to match 4-product global mean
gcb4prod=1.136
lulucnet=lulucnet*gcb4prod/lulucglb
luluctrp=luluctrp*gcb4prod/lulucglb
lulucset=lulucset*gcb4prod/lulucglb
lulucglb=gcb4prod

#> c(lulucset,luluctrp,lulucnet)
#[1] 0.07216015 0.90243476 0.16140509

## Friedlingstein et al., Nature, 2026 reported an increase in Eluc of +0.11 +/- 0.04 for 2014-2023 ("delta-L")
## adjust fluxes proportionally to match new global mean
friedglb=lulucglb+0.11
lulucnet=lulucnet*friedglb/lulucglb
luluctrp=luluctrp*friedglb/lulucglb
lulucset=lulucset*friedglb/lulucglb
lulucglb=friedglb

## reduction by 0.12 and increase by 0.11 should has almost no net effect
#> c(lulucset,luluctrp,lulucnet)
#[1] 0.07914749 0.98981841 0.17703410

## for uncertainty use sqrt(flux fraction) from GCB2024
## GCB2024 uses expert judgement to assign 0.7 PgC/yr uncertainty globally
# since all fluxes are positive, do not need to use abs values
unclulucglb=0.7
## add in 0.04 from Friedlingstein adjustment in quadrature
unclulucglb=(unclulucglb^2+0.04^2)^0.5
#unclulucglb
#[1] 0.7011419
## almost no effect

print('LULUC')
print('(BLUE and LUCE divided at 20S/20N for 2016-2018, scaled to match 4-product global mean 2016-2018 and adjusted by Friedlingstein et al. (Nature, 2026) for 2014-2023)')
print(paste('GLB:',round(lulucglb,2),'+/-',round(unclulucglb,2)))
unclulucnet=unclulucglb*(lulucnet/lulucglb)^0.5
print(paste('NET:',round(lulucnet,2),'+/-',round(unclulucnet,2)))
uncluluctrp=unclulucglb*(luluctrp/lulucglb)^0.5
print(paste('TRP:',round(luluctrp,2),'+/-',round(uncluluctrp,2)))
unclulucset=unclulucglb*(lulucset/lulucglb)^0.5
print(paste('SET:',round(lulucset,2),'+/-',round(unclulucset,2)))
print((unclulucnet^2+uncluluctrp^2+unclulucset^2)^0.5)


## Rivers:

# GCB2024: 
# "The global fCO2-based flux estimates were adjusted to remove the pre-industrial ocean source of CO2 to the atmosphere of 0.65 ± 0.3 GtC yr−1 from river input to the ocean (Regnier et al., 2022) in order to satisfy our definition of SOCEAN (Hauck et al., 2020). The river flux adjustment was distributed over the latitudinal bands using the regional distribution of Lacroix et al. (2020; north: 0.14 GtC yr−1 ; tropics: 0.42 GtC yr−1 ; south: 0.09 GtC yr−1 ).
## However, to correct DGVM fluxes, we want the distribution of the land uptake, not the ocean outgassing
### these are available in the GCB 1x1 gridded inversion file: GCP2024_inversions_1x1_version1_2_20241021.nc
# double river_adjustment_land_scaled[longitude,latitude]   (Chunking: [360,180])  (Compression: shuffle,level 4)
#     long_name: river adjustment for the land but scaled to match the numbers provided for the latitudes: 0.30 N, 0.27
# from netcdf header:
# Lateral river flux adjustment on land as provided by Ronny Lauerwald (The file is based on GlobalNEWS2 for organic C and the weathering CO2 sink after Hartmann et al. 2009 as used in Zscheischler et al 2017. But in this version, the organic C loads after GlobalNEWS are twice rescaled: 1) to the latitudinal pattern from Resplandy et al. (2018 NatGeo) and 2) to a synthesis of global estimates of organic C exports of about 500 Tg C/yr (for this you could for the time being cite Regnier et al. 2013, Nat Geo).).
### processed by GCB2024/flux_proc_gcb2024.r and ouput as zonal files, read back in here

library('abind')

## Options:
gcbfluxdir='LOCALDATA'
gcbfluxfile='GCP2024_inversions_1x1_version1_2_20241021.nc'

cacalclist=c( # slat, nlat, limtxt
20,90,'20N-90N',
-20,20,'20S-20N',
-90,-20,'90S-20S'
)

caspecs=matrix(cacalclist,byrow=T,ncol=3)

# open inversion flux file
fluxnc=nc_open(paste(gcbfluxdir,'/',gcbfluxfile,sep=''))
# double river_adjustment_land[longitude,latitude]   (Chunking: [360,180])  (Compression: shuffle,level 4)
#     long_name: river adjustment for the land, see summary in the netcdf header for more information

area=ncvar_get(fluxnc,'cell_area')
latarea=apply(area,2,mean)
flat=ncvar_get(fluxnc,'latitude')

riverfluxes_lonlat=ncvar_get(fluxnc,'river_adjustment_land_scaled') # constant and same for all models (lon x lat)
riverfluxes_lonlat=sweep(riverfluxes_lonlat, MARGIN = 2, STATS = latarea, FUN = "*")

# trimmed by lat on each loop so save to reset
riverfluxsave=riverfluxes_lonlat

for(i in c(1:nrow(caspecs))){

        riverflux=riverfluxsave
        print(caspecs[i,])
        latlim=as.numeric(caspecs[i,1:2])
        limtxt=caspecs[i,3]

        # mask for region
        if(grepl('X',limtxt)){ # lump together SH and NH
                riverflux[,abs(flat)<latlim[1]|abs(flat)>=latlim[2]]=0
        } else {
                riverflux[,flat<latlim[1]|flat>=latlim[2]]=0
        }

        # sum for each month
        regriverflux=sum(riverflux,na.rm=T)
        # write out river flux
        write(regriverflux,paste('GCB2024/riverflux_',limtxt,'.txt',sep=''))

} # loop on caspecs

# read back in:
rivnet=scan('GCB2024/riverflux_20N-90N.txt')*-1
rivtrp=scan('GCB2024/riverflux_20S-20N.txt')*-1
rivset=scan('GCB2024/riverflux_90S-20S.txt')*-1
#print(round(c(rivset+rivtrp+rivnet,rivset,rivtrp,rivnet),2))
#[1] -0.65 -0.09 -0.21 -0.35

rivglb=rivset+rivtrp+rivnet

## for uncertainty use sqrt(flux fraction) from GCB2024
## since all fluxes are positive, do not need to use abs values
uncrivglb=0.3

print('Rivers')
print('(GCB2024 global total of 0.65 scaled by land uptake)')
print(paste('GLB:',round(rivglb,2),'+/-',round(uncrivglb,2)))
uncrivnet=uncrivglb*(rivnet/rivglb)^0.5
print(paste('NET:',round(rivnet,2),'+/-',round(uncrivnet,2)))
uncrivtrp=uncrivglb*(rivtrp/rivglb)^0.5
print(paste('TRP:',round(rivtrp,2),'+/-',round(uncrivtrp,2)))
uncrivset=uncrivglb*(rivset/rivglb)^0.5
print(paste('SET:',round(rivset,2),'+/-',round(uncrivset,2)))
print((uncrivnet^2+uncrivtrp^2+uncrivset^2)^0.5)
uncrivsettrp=uncrivglb*((rivset+rivtrp)/rivglb)^0.5
print(paste('SET+T:',round(rivset+rivtrp,2),'+/-',round(uncrivsettrp,2)))


# combine
adjglb=nbpglb+rssglb+lulucglb+lceglb+rivglb
adjnet=nbpnet+rssnet+lulucnet+lcenet+rivnet
adjtrp=nbptrp+rsstrp+luluctrp+lcetrp+rivtrp
adjset=nbpset+rssset+lulucset+lceset+rivset
uncadjglb=(uncnbpglb^2+uncrssglb^2+unclulucglb^2+unclceglb^2+uncrivglb^2)^0.5
uncadjnet=(uncnbpnet^2+uncrssnet^2+unclulucnet^2+unclcenet^2+uncrivnet^2)^0.5
uncadjtrp=(uncnbptrp^2+uncrsstrp^2+uncluluctrp^2+unclcetrp^2+uncrivtrp^2)^0.5
uncadjset=(uncnbpset^2+uncrssset^2+unclulucset^2+unclceset^2+uncrivset^2)^0.5


# write out results
print('Table S5 (DGVM adjustements)')
print(paste('GCB2024 DGVMs$^1$ \\newline (June 2016--May 2018) & ',round(nbpset,2),' $\\pm$ ',round(uncnbpset,2),' & ',round(nbptrp,2),' $\\pm$ ',round(uncnbptrp,2),' & ',round(nbpnet,2),' $\\pm$ ',round(uncnbpnet,2),' \\',sep=''))
print(paste('RSS$^2$ \\newline (2015--2024) & ',round(rssset,2),' $\\pm$ ',round(uncrssset,2),' & ',round(rsstrp,2),' $\\pm$ ',round(uncrsstrp,2),' & ',round(rssnet,2),' $\\pm$ ',round(uncrssnet,2),' \\',sep=''))
print(paste('LULUC$^3$ \\newline (2016--2018) & ',round(lulucset,2),' $\\pm$ ',round(unclulucset,2),' & ',round(luluctrp,2),' $\\pm$ ',round(uncluluctrp,2),' & ',round(lulucnet,2),' $\\pm$ ',round(unclulucnet,2),' \\',sep=''))
print(paste('LCE$^4$ \\newline (2014--2023) & ',round(lceset,2),' $\\pm$ ',round(unclceset,2),' & ',round(lcetrp,2),' $\\pm$ ',round(unclcetrp,2),' & ',round(lcenet,2),' $\\pm$ ',round(unclcenet,2),' \\',sep=''))
print(paste('Rivers$^5$ \\newline (pre-industrial) & ',round(rivset,2),' $\\pm$ ',round(uncrivset,2),' & ',round(rivtrp,2),' $\\pm$ ',round(uncrivtrp,2),' & ',round(rivnet,2),' $\\pm$ ',round(uncrivnet,2),' \\',sep=''))
print(paste('Adjusted land-to-air flux & ',round(adjset,2),' $\\pm$ ',round(uncadjset,2),' & ',round(adjtrp,2),' $\\pm$ ',round(uncadjtrp,2),' & ',round(adjnet,2),' $\\pm$ ',round(uncadjnet,2),' \\',sep=''))
print('$^1$Mean of 20 models from Friedlingstein et al. 2025') 
print('$^2$Replaced Sources and Sinks. Mean of 7 models from O\'Sullivan et al. 2026') 
print('$^3$Mean of BLUE and LUCE adjusted to match GCB2024 4-product global mean and for transient carbon densities from Friedlingstein et al. 2026') 
print('$^4$Lateral Carbon Export. From Friedlingstein et al. 2026 and Tian et al. 2023')
print('$^5$From Friedlingstein et al. 2025 and Regnier et al. 2013 scaled by uptake on land')
print('remove slashes from \\newline and \\pm, and add zeros')


print(round(c(adjglb,adjset,adjtrp,adjnet),2))
