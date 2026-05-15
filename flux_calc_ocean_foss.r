## Program to read in and process independent air-sea CO2 and Fossil-fuel emission estimates


library('ncdf4')

# Ocean fluxes from pCO2 data products in GCB2024, downloaded from:
#Hauck, J., Mayot, N., Landschützer, P., & Jersild, A. (2025). Global Carbon Budget 2024, surface ocean fugacity of CO2 (fCO2) and air-sea CO2 flux of individual global ocean biogeochemical models and surface ocean fCO2-based data-products [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14639761
# version2 
pco2dir='LOCALDATA/GCBOCEAN'
files=list.files(pco2dir,pattern='_dataprod')
models=unlist(lapply(strsplit(files,'_'),function(x) x[3]))

# not including UoEX-UEPFFNU
# "A ninth fCO2 product (UExP-FFN-U) is shown but is not included in the ensemble average as it differs from the other products by adjusting the flux to a cool, salty ocean surface skin."
files=files[models!='UoEX-UEPFFNU']
models=models[models!='UoEX-UEPFFNU']

# GCB 2024 scales CMEMS-LSCE-FFNN, NIES-ML3, LDEO-HPD, CSIR-ML6 for coastal regions (note SM has a typo, LDEO was scaled not UoEX)
toscale=c('CMEMS-LSCE-FFNN','NIES-ML3','LDEO-HPD','CSIR-ML6')
toscaleshortnames=c('CMEMS','NIES','LDEO','CSIR')
## CMEMS-LSCE-FFNN appears to end in 2009
# Peter L. had monthly 'true' ice free areas, MaxArea.mat = monthly values for entire timeseries (396 months 1993-2021 inclusive) by 4 regions (glob, north, tropics, south) dividing at 30 degrees
# Peter L. had annual scale factors, ScalingFactors.mat = annual scale factors for 8 products and 4 regions (glob, north, tropics, south) dividing at 30 degrees
## missing JENA-MLS , but Jena is one of the products he says does not need scaling
# for the 4 that do:
# $ CMEMSPerc: num [1, 1:4] 0.944 0.885 0.941 0.983
# $ CSIRPerc : num [1, 1:4] 0.962 0.856 0.976 0.993
# $ LDEOPerc : num [1, 1:4] 0.935 0.867 0.932 0.98
# $ NIESPerc : num [1, 1:4] 0.949 0.894 0.947 0.985
# and 4 of the 5 that don't:
# $ ETHPerc  : num [1, 1:4] 1 1.01 1 1
# $ JMAPerc  : num [1, 1:4] 1.01 1.06 1.01 1.01
# $ UoEXPerc : num [1, 1:4] 1.008 1.043 0.999 1.005
# $ VLIZPerc : num [1, 1:4] 0.991 0.98 0.991 0.999
# Requested and Peter added monthly product areas (denominator in equation S4), Denominator_*product*.mat = monthly values for entire timeseries (396 months 1993-2021 inclusive) by 4 regions (glob, north, tropics, south) dividing at 30 degrees
## product areas are smaller than 'true' areas so scalefactors being <1 means (corrected flux) = (product flux)/scalefactor = (true area) / (product area) * (product flux)
# file names: Denominator_CMEMS.mat  Denominator_CSIR.mat  Denominator_Jena.mat  Denominator_JMA.mat  Denominator_LDEO.mat  Denominator_NIES.mat  Denominator_OceanSODA.mat  Denominator_SOMFFN.mat  Denominator_UoEX.mat
## So, to adjust these 4 products, 
#1) read in MaxArea.mat
#2) select variable North1, South1, or Tropics1 and select June 2016-May 2018
#3) read in Denominator_*product*.mat for that product
#4) select variable area.north, area.south, or area.tropics and select June 2016-May 2018
#4) divide maxarea by denominator and multiply by flux
## to do 20 degree, need gridded monthly maxareas and gridded monthly product areas
library(R.matlab)
maxareafile=paste(pco2dir,'/dataprod_areas/MaxArea.mat',sep='')
maxarea=readMat(maxareafile)
areadt=seq.POSIXt(ISOdatetime(1991,1,15,0,0,0,tz='UTC'),ISOdatetime(2023,12,15,0,0,0,tz='UTC'),by='month')
areazacsel=difftime(areadt,ISOdatetime(2016,6,1,0,0,0,tz='UTC'))>0&difftime(areadt,ISOdatetime(2018,6,1,0,0,0,tz='UTC'))<0

## reg variables are for 30S and 30N
# *regall uses provided 4 (30 deg divided) regions
glbregall=NULL
netregall=NULL
trpregall=NULL
setregall=NULL
soregall=NULL
emtregall=NULL # ET-T
# *all uses gridded to get 20 deg divisions
glball=NULL
netall=NULL
trpall=NULL
setall=NULL
soall=NULL
emtall=NULL
latdiv=c(-20,20)
areaall=NULL
for(i in c(1:length(files))){
	print(files[i])
	nc=nc_open(paste(pco2dir,'/',files[i],sep=''))
	fgco2=ncvar_get(nc,'fgco2')
	fgco2_reg=ncvar_get(nc,'fgco2_reg')
	print(dim(fgco2_reg))
	lat=ncvar_get(nc,'lat')
	area=ncvar_get(nc,'area')
	# SODA area is 16X too small
	if(models[i]=="OceanSODA-ETHZv2") area=area*16
	areaall=c(areaall,sum(area,na.rm=T))
	# double area[lon,lat]   (Contiguous storage) # units: m2
	# float fgco2[lon,lat,time]   (Contiguous storage) # units: mol/m2/s
	time=ncvar_get(nc,'time')
        # units: days since 1959-01-01
	dt=ISOdatetime(1959,1,1,0,0,0,tz='UTC')+(time+15)*86400
	print(range(dt))
	zacsel=difftime(dt,ISOdatetime(2016,6,1,0,0,0,tz='UTC'))>0&difftime(dt,ISOdatetime(2018,6,1,0,0,0,tz='UTC'))<0
#	zacsel=difftime(dt,ISOdatetime(2023,1,1,0,0,0,tz='UTC'))>0&difftime(dt,ISOdatetime(2024,1,1,0,0,0,tz='UTC'))<0
#	zacsel=difftime(dt,ISOdatetime(2016,1,1,0,0,0,tz='UTC'))>0&difftime(dt,ISOdatetime(2017,1,1,0,0,0,tz='UTC'))<0
#	zacsel=difftime(dt,ISOdatetime(2000,1,1,0,0,0,tz='UTC'))>0&difftime(dt,ISOdatetime(2020,1,1,0,0,0,tz='UTC'))<0 # for comparing to Randerson et al
	glb=NULL
	net=NULL
	trp=NULL
	set=NULL
	so=NULL
	emt=NULL
	for(j in c(1:length(time))[zacsel]){
		fgco2tot=fgco2[,,j]*area
		glb=c(glb,sum(fgco2tot,na.rm=T)*12/1E15*86400*365) # not accounting for leap year (2016)
		net=c(net,sum(fgco2tot[,lat>latdiv[2]],na.rm=T)*12/1E15*86400*365)
		trp=c(trp,sum(fgco2tot[,lat<latdiv[2]&lat>latdiv[1]],na.rm=T)*12/1E15*86400*365)
		set=c(set,sum(fgco2tot[,lat<latdiv[1]],na.rm=T)*12/1E15*86400*365)
		so=c(so,sum(fgco2tot[,lat<(-45)],na.rm=T)*12/1E15*86400*365)
	}
	## adjust for coastal areas
	if(is.element(models[i],toscale)){
		file=paste(pco2dir,'/dataprod_areas/Denominator_',toscaleshortnames[which(toscale==models[i])],'.mat',sep='')
		denom=readMat(file)
		glb=maxarea$Global1[areazacsel]/denom$area.global[areazacsel]*glb
		net=maxarea$North1[areazacsel]/denom$area.north[areazacsel]*net
		trp=maxarea$Tropics[areazacsel]/denom$area.tropics[areazacsel]*trp
		set=maxarea$South1[areazacsel]/denom$area.south[areazacsel]*set
	}
	emt=net+set-trp
	if(models[i]=='JENA-MLS'){ # JENA region fluxes are annual (OK, because not using these)
		glbregall=c(glbregall,NA) # can't do ATom period with whole years
		netregall=c(netregall,NA)
		trpregall=c(trpregall,NA)
		setregall=c(setregall,NA)
		soregall=c(soregall,NA)
		emtregall=c(emtregall,NA)
	} else if(which(dim(fgco2_reg)==4)==2){ # some are rotated
		glbregall=c(glbregall,mean(fgco2_reg[zacsel,1],na.rm=T))
		netregall=c(netregall,mean(fgco2_reg[zacsel,2],na.rm=T))
		trpregall=c(trpregall,mean(fgco2_reg[zacsel,3],na.rm=T))
		setregall=c(setregall,mean(fgco2_reg[zacsel,4],na.rm=T))
		soregall=c(soregall,mean(fgco2_reg[zacsel,4],na.rm=T))
		emtregall=c(emtregall,mean(fgco2_reg[zacsel,2]+fgco2_reg[zacsel,4],na.rm=T)-fgco2_reg[zacsel,3])
	} else {
		glbregall=c(glbregall,mean(fgco2_reg[1,zacsel],na.rm=T))
		netregall=c(netregall,mean(fgco2_reg[2,zacsel],na.rm=T))
		trpregall=c(trpregall,mean(fgco2_reg[3,zacsel],na.rm=T))
		setregall=c(setregall,mean(fgco2_reg[4,zacsel],na.rm=T))
		soregall=c(soregall,mean(fgco2_reg[4,zacsel],na.rm=T))
		emtregall=c(emtregall,mean(fgco2_reg[2,zacsel]+fgco2_reg[4,zacsel],na.rm=T)-fgco2_reg[3,zacsel])
	}
	glball=c(glball,mean(glb,na.rm=T))
	netall=c(netall,mean(net,na.rm=T))
	trpall=c(trpall,mean(trp,na.rm=T))
	setall=c(setall,mean(set,na.rm=T))
	soall=c(soall,mean(so,na.rm=T))
	emtall=c(emtall,mean(emt,na.rm=T))
	nc_close(nc)
}
print(round(glball,2))
print(round(netall,2))
print(round(trpall,2))
print(round(setall,2))
print(round(soall,2))
print(round(emtall,2))
mean(glball)
mean(netall)
mean(trpall)
mean(setall)
mean(soall)
mean(emtall)
#[1] 2.422739
#[1] 1.338026
#[1] -0.564582
#[1] 1.649295
sd(glball)
sd(netall)
sd(trpall)
sd(setall)
sd(soall)
sd(emtall)
#[1] 0.2962539
#[1] 0.09851352
#[1] 0.1249872
#[1] 0.130024
glboce=-1*mean(glball)
netoce=-1*mean(netall)
trpoce=-1*mean(trpall)
setoce=-1*mean(setall)
sooce=-1*mean(soall)
emtoce=-1*mean(emtall)

## to add in systemmatic errors from GCB 2024, ignore river term (accounted elsewhere) and treat the rest as if they scale by flux for regions (100% correlated)
# 0.2 for obs
# 0.2 for gas exchange
# 0.1 for winds
# 0.2 for mapping

uncglboce=(sd(glball)^2+0.2^2+0.2^2+0.1^2+0.2^2)^0.5
print(paste('GLB:',round(glboce,2),'+/-',round(uncglboce,2)))

frat=(abs(netoce)/(abs(netoce)+abs(trpoce)+abs(setoce)))^0.5 # div by sum of abs fluxes rather than glb
uncnetoce=(sd(netall)^2+(0.2*frat)^2+(0.2*frat)^2+(0.1*frat)^2+(0.2*frat)^2)^0.5
print(paste('NET:',round(netoce,2),'+/-',round(uncnetoce,2)))

frat=(abs(trpoce)/(abs(netoce)+abs(trpoce)+abs(setoce)))^0.5 # div by sum of abs fluxes rather than glb
unctrpoce=(sd(trpall)^2+(0.2*frat)^2+(0.2*frat)^2+(0.1*frat)^2+(0.2*frat)^2)^0.5
print(paste('TRP:',round(trpoce,2),'+/-',round(unctrpoce,2)))

frat=(abs(setoce)/(abs(netoce)+abs(trpoce)+abs(setoce)))^0.5 # div by sum of abs fluxes rather than glb
uncsetoce=(sd(setall)^2+(0.2*frat)^2+(0.2*frat)^2+(0.1*frat)^2+(0.2*frat)^2)^0.5
print(paste('SET:',round(setoce,2),'+/-',round(uncsetoce,2)))

sd(netall+trpall+setall)

# using 30S/N and 2023 got
#[1] 2.351515
#[1] 1.148086
#[1] -0.08486621
#[1] 1.288295
# adding in rivers (Globe: 0.63, North: 0.14 GtC yr-1, Tropics: 0.42 GtC yr-1, South: 0.09 GtC yr-1, see GCB 2024 paper)
# gives GCB Fig. 11 (3) and 14 (1.3, .3, 1.4)



# Fossil-fuel emissions from GridFED, Jones et al., 2024 (update to 2021)
## for 2023, got 10.04 for globe vs. 10.1 reported in GCB 2024
fossdir='LOCALDATA/GCBFOSS'
glbff=0
netff=0
trpff=0
setff=0
emtff=0
allsel=0
for(year in c(2016:2018)){
	nc=nc_open(paste(fossdir,'/GCP-GridFEDv2024.0_',year,'.nc',sep=''))
	#        float CO2/TOTAL[lon,lat,time]   (Chunking: [1200,600,1])
	#            units: kg CO2 cell-1 month-1
	ff=ncvar_get(nc,'CO2/TOTAL')
	lat=ncvar_get(nc,'lat')
	time=ncvar_get(nc,'time')
	dt=ISOdatetime(year,1,1,0,0,0,tz='UTC')+(time+15)*86400
	zacsel=difftime(dt,ISOdatetime(2016,6,1,0,0,0,tz='UTC'))>0&difftime(dt,ISOdatetime(2018,6,1,0,0,0,tz='UTC'))<0
	glbff=glbff+sum(ff[,,zacsel])
	netff=netff+sum(ff[,lat>20,zacsel])
	trpff=trpff+sum(ff[,lat>(-20)&lat<20,zacsel])
	setff=setff+sum(ff[,lat<(-20),zacsel])
	emtff=emtff+sum(ff[,lat<(-20),zacsel])+sum(ff[,lat>20,zacsel])-sum(ff[,lat>(-20)&lat<20,zacsel])
	allsel=allsel+sum(zacsel)
	nc_close(nc)
}
glbff=glbff/1E12*12/44*12/allsel # to PgC and /year
netff=netff/1E12*12/44*12/allsel # to PgC and /year
trpff=trpff/1E12*12/44*12/allsel # to PgC and /year
setff=setff/1E12*12/44*12/allsel # to PgC and /year
emtff=emtff/1E12*12/44*12/allsel # to PgC and /year
print(glbff)
print(netff)
print(trpff)
print(setff)
print(emtff)
#[1] 9.571133
#[1] 8.187921
#[1] 0.9758878
#[1] 0.4073244

# for uncertainty, follow GCB 2024 and Jones et al. 2021
# assume 5% for NET except for 10% on China
# assume 10% for TRP
# assume 10% for SET except for Australia and New Zealand
uncglbff=0.05*glbff
print(paste('GLB:',round(glbff,2),'+/-',round(uncglbff,2)))
china=2.73
uncnetff=0.05*(netff-china)+0.1*china
uncnetpercent=uncnetff/netff
print(paste('NET:',round(netff,2),'+/-',round(uncnetff,2),'(',round(uncnetpercent*100,2),'%)'))
unctrpff=0.1*trpff
unctrppercent=unctrpff/trpff
print(paste('TRP:',round(trpff,2),'+/-',round(unctrpff,2),'(',round(unctrppercent*100,2),'%)'))
aust=0.11; nz=0.01
uncsetff=0.1*(setff-aust-nz)+0.05*(aust+nz)
uncsetpercent=uncsetff/setff
print(paste('SET:',round(setff,2),'+/-',round(uncsetff,2),'(',round(uncsetpercent*100,2),'%)'))

uncemtff=(uncnetff^2+uncsetff^2+unctrpff^2)^0.5 ## assume uncertainties are uncorrelated bw net, trp, and set
uncemtpercent=uncemtff/emtff
print(paste('EMT:',round(emtff,2),'+/-',round(uncemtff,2),'(',round(uncemtpercent*100,2),'%)'))

# from https://globalcarbonbudgetdata.org/downloads/jGJH0-data/National_Fossil_Carbon_Emissions_2024v1.0.xlsx
# 0.11276934 # =(7/12*L179+L180+5/12*L181)/2/1000
# 2.729338268
# 0.009617744 

#[1] "GLB: 9.57 +/- 0.48"
#[1] "NET: 8.19 +/- 0.55"
#[1] "TRP: 0.98 +/- 0.1"
#[1] "SET: 0.41 +/- 0.03"

#ATom / v10 MIP ECFC Total Flux
#-1.82 ± 0.22 1.92 ± 0.39 4.60 ± 0.28
#Fossil-fuel Emission (GridFED)
#0.41 ± 0.03 0.98 ± 0.10 8.19 ± 0.55
#Ocean Flux (Hauck et al., 2025)
#-1.65 ± 0.21 0.56 ± 0.14 -1.34 ± 0.17
print(paste('NET:',4.60-8.19+1.34,'+/-',round((0.28^2+0.55^2+0.17^2)^0.5,2)))
print(paste('TRP:',1.92-0.98-0.56,'+/-',round((0.37^2+0.10^2+0.14^2)^0.5,2)))
print(paste('SET:',-1.81-0.41+1.65,'+/-',round((0.22^2+0.03^2+0.21^2)^0.5,2)))
#[1] "NET: -2.25 +/- 0.64"
#[1] "TRP: 0.38 +/- 0.41"
#[1] "SET: -0.57 +/- 0.31"


print('Table 1 (results)')
## first part from flux_conc_ecfc.r
print(paste('pCO$_2$-based \\newline sea-to-air flux \\cite{Hauck2025} & ',round(setoce,2),' $\\pm$ ',round(uncsetoce,2),' & ',round(trpoce,2),' $\\pm$ ',round(unctrpoce,2),' & ',round(netoce,2),' $\\pm$ ',round(uncnetoce,2),' \\',sep=''))
print(paste('Canonical fossil-fuel \\newline emission \\cite{Jones2024} & ',round(setff,2),' $\\pm$ ',round(uncsetff,2),' & ',round(trpff,2),' $\\pm$ ',round(unctrpff,2),' & ',round(netff,2),' $\\pm$ ',round(uncnetff,2),' \\',sep=''))


## from flux_conc_ecfc.r
#netecfc=4.594; uncnetecfc=0.282
#trpecfc=1.925; unctrpecfc=0.391
#setecfc=-1.821; uncsetecfc=0.232
#netecfc=4.588; uncnetecfc=0.283
#trpecfc=1.922; unctrpecfc=0.387
#setecfc=-1.819; uncsetecfc=0.234
netecfc=4.592; uncnetecfc=0.284
trpecfc=1.919; unctrpecfc=0.385
setecfc=-1.82; uncsetecfc=0.23


netres=netecfc-netff-netoce
trpres=trpecfc-trpff-trpoce
setres=setecfc-setff-setoce
uncnetres=(uncnetecfc^2+uncnetff^2+uncnetoce^2)^0.5
unctrpres=(unctrpecfc^2+unctrpff^2+unctrpoce^2)^0.5
uncsetres=(uncsetecfc^2+uncsetff^2+uncsetoce^2)^0.5
print(paste('\textbf{Canonical implied \\newline land-toair flux} & \textbf{',round(setres,2),' $\\pm$ ',round(uncsetres,2),'} & \textbf{',round(trpres,2),' $\\pm$ ',round(unctrpres,2),'} & \textbf{',round(netres,2),' $\\pm$ ',round(uncnetres,2),'} \\',sep=''))
