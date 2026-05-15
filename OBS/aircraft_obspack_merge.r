# Program to merge variables camp, prof, pressure, and theta from aircraft campaign merge files onto NOAA GlobalView+ Obspack IDs
# Outputs e.g. OBS/ATOM_obspack_GLOBALVIEWplus_v6.1_merge.txt, which is then used by V10MIP/model_aircraft_merge.r 

## note on scales: all 10-sec merge product CO2 data read in here is reported on the WMO X2007 scale, CO2 data from ObsPack files is not read in here

## To do:
# get Mthetae from official repo or remove option

aomerge<-function(gvpvernum='6.1'){

	library('ncdf4')

	## Set ObsPack Globalview+ name
	opname=paste('GLOBALVIEWplus_v',gvpvernum,sep='')

	## Read in N2O data for stratospheric filter reference:
	glbn2ofile=url('ftp://aftp.cmdl.noaa.gov/products/trends/n2o/n2o_annmean_gl.txt')
	temp=readLines(glbn2ofile)
	hlines<-which(temp=='# year     mean      unc')
	glbn2o=read.table(glbn2ofile,skip=hlines,header=F,stringsAsFactors=F)
	colnames(glbn2o)=c('year','n2o','unc')

	## Set strat flag cutoffs based upon Jin et al.: https://acp.copernicus.org/preprints/acp-2020-841/acp-2020-841.pdf
	stratcoh2o=50
	stratcoo3=150
	stratcon2o=319

	## Read in ATom file, calculate stratosphere flag, select latitudes and columns, and subset
	atommergedir='../LOCALDATA/ATOM' # these are version 2.0 (21-08-26) from https://doi.org/10.3334/ORNLDAAC/1925
	atommergefiles=c('MER10_DC8_ATom-1.nc','MER10_DC8_ATom-2.nc','MER10_DC8_ATom-3.nc','MER10_DC8_ATom-4.nc')

	atomvar=c('UTC_Start','Flight_Date','time','DLH-H2O/H2O_DLH','UCATS-H2O/H2O_UWV','QCLS-CH4-CO-N2O/N2O_QCLS','GCECD/N2O_PECD','UCATS-GC/N2O_UCATS','NOyO3-O3/O3_CL','UCATS-O3/O3_UCATS','MMS/G_ALT','RF','prof.no','MMS/P','MMS/POT','MMS/G_LAT','MMS/G_LONG','MMS/HDG','NOAA-Picarro/CO2_NOAA','QCLS-CO2/CO2_QCLS','AO2/CO2_AO2','CO2.X','NOAA-Picarro/CH4_NOAA','QCLS-CH4-CO-N2O/CH4_QCLS','UCATS-GC/CH4_UCATS','GCECD/CH4_PECD','GCECD/SF6_PECD','UCATS-GC/SF6_UCATS')

	atommerge=NULL
	for(i in c(1:4)){
		atomnc=nc_open(paste(atommergedir,'/',atommergefiles[i],sep=''))
		count=length(ncvar_get(atomnc,'UTC_Start'))
		campdata=NULL
		for(var in atomvar){
			if(i==1&var=='UCATS-H2O/H2O_UWV'){ # no UCATS H2O on ATom-1
				campdata=cbind(campdata,rep(NA,count))
			} else {
				campdata=cbind(campdata,ncvar_get(atomnc,var))
			}
		}
		campdata=cbind(campdata,rep(i,count)) # A.no
		nc_close(atomnc)
		atommerge=rbind(atommerge,campdata)
	}
	atommerge=data.frame(atommerge,stringsAsFactors=F)
	names(atommerge)=c(gsub('.*/','',atomvar),'A.no')

	# add time variables
	atommerge$Year=as.numeric(substr(atommerge$Flight_Date,1,4))
	atommerge$Month=as.numeric(substr(atommerge$Flight_Date,5,6))
	atommerge$Day=as.numeric(substr(atommerge$Flight_Date,7,8)) # Flight_Date is day of takeoff
	atommerge$UTC_Mean=atommerge$UTC_Start+5
	atomdt=as.POSIXlt(ISOdatetime(atommerge$Year,atommerge$Month,atommerge$Day,0,0,0,tz='UTC')+atommerge$UTC_Mean,tz='UTC') # use UTC_Mean to match actual sample time
	atommerge$DOY=as.POSIXlt(atomdt)$yday+1 # atomdt corresponds to UTC_Mean
	atommerge$DOY[atommerge$UTC_Mean>86400]=atommerge$DOY[atommerge$UTC_Mean>86400]-1 ## do not want DOY to jump a whole day over midnight
	atommerge$jday2016=atommerge$time/86400
	atommerge$YYYYMMDD=strftime(as.POSIXlt(atomdt),format='%Y%m%d',tz='UTC')

	atomdt=as.POSIXlt(ISOdatetime(atommerge$Year,atommerge$Month,atommerge$Day,0,0,0,tz='UTC')+atommerge$UTC_Start,tz='UTC') # now use UTC_Start to merge with ObsPack
	atommerge$Year=as.POSIXlt(atomdt)$year+1900 # reset to account for flights over midnight UTC
	atommerge$Month=as.POSIXlt(atomdt)$mon+1
	atommerge$Day=as.POSIXlt(atomdt)$mday
	atommerge$Hour=as.POSIXlt(atomdt)$hour
	atommerge$Min=as.POSIXlt(atomdt)$min
	atommerge$Sec=as.POSIXlt(atomdt)$sec

	# add strat flag
	atommerge$strat=rep(0,nrow(atommerge)) # 0 means trop
	h2oref=atommerge$H2O_DLH; h2oref[is.na(h2oref)]=atommerge$H2O_UWV[is.na(h2oref)]; h2oref[is.na(h2oref)]=0 # if H2O missing treat as if potentially strat
	n2oref=atommerge$N2O_QCLS; n2oref[is.na(n2oref)]=atommerge$N2O_PECD[is.na(n2oref)]; n2oref[is.na(n2oref)]=atommerge$N2O_UCATS[is.na(n2oref)]
	n2oref=n2oref-(approx(glbn2o$year+0.5,glbn2o$n2o,atommerge$Year+atommerge$DOY/365)$y-glbn2o$n2o[glbn2o$year==2009])
	n2oref[is.na(n2oref)]=400 # if N2O missing do not use for filter
	o3ref=atommerge$O3_CL; o3ref[is.na(o3ref)]=atommerge$O3_UCATS[is.na(o3ref)]; o3ref[is.na(o3ref)]=0 # if O3 missing do not use for filter
	atommerge$strat[h2oref<stratcoh2o&(o3ref>stratcoo3|n2oref<stratcon2o|(o3ref==0&n2oref==400))]=1 # if either o3 or n2o criteria are met, or if both are missing, consider strat
	atommerge$strat[h2oref==0&o3ref==0&n2oref==400&atommerge$G_ALT<8000]=0 # if all 3 missing assume < 8 km is trop
	colsel=c('Year','Month','Day','Hour','Min','Sec','A.no','RF','prof.no','P','POT','strat','HDG','m_theta_e','G_LAT','G_LONG','G_ALT','CO2_QCLS','CO2_NOAA','CO2_AO2','CO2.X')
	atommerge=atommerge[,is.element(colnames(atommerge),colsel)]
	atommerge=atommerge[,match(colsel,names(atommerge),nomatch=F)] # reorder
	names(atommerge)=c('year','mon','day','hour','min','sec','camp','flt','prof','pressure','theta','strat','hdg','lat','lon','alt','CO2_QCLS','CO2_NOAA','CO2_AO2','CO2.X')
	write(names(atommerge),ncol=ncol(atommerge),'ATom_trimmed.txt')
	write(t(atommerge),'ATom_trimmed.txt',ncol=ncol(atommerge),append=T)

	## After writing out version with CO2 values, trim back further to exclude these
	atommerge=atommerge[,!grepl('CO2',names(atommerge))]

	## Merge obspack_id and obspack_num onto 10-sec merge file data by c('year','mon','day','hour','min','sec')
	camp='ATom'
	ncfile='co2_tom_aircraft-insitu_1_allvalid.nc'
	print(paste(opname,camp))
	## Read in observed obspack files
	ncin=nc_open(paste('../LOCALDATA/',ncfile,sep=''))
	ncdat=data.frame(cbind(t(ncvar_get(ncin,'time_components')),ncvar_get(ncin,'obspack_id'),ncvar_get(ncin,'obspack_num'))) ; colnames(ncdat)=c('year','mon','day','hour','min','sec','obspack_id','obspack_num')
	if(as.numeric(gvpvernum)>=7) ncdat$sec=as.numeric(ncdat$sec)-5 # GV+ 7.0 and after use mid instead of start time from 10-sec merge
	mergefile=get(paste(tolower(camp),'merge',sep=''))
	## number of data points greater in merge files because of test flights, missing data, calibrations, etc.
	mrgdat=merge(mergefile,ncdat,by=c('year','mon','day','hour','min','sec'))
	print(dim(mergefile))
	print(dim(ncdat))
	print(dim(mrgdat))
	write(colnames(mrgdat),paste(camp,'_obspack_',opname,'_merge.txt',sep=''),ncol=ncol(mrgdat))
	write(t(mrgdat),paste(camp,'_obspack_',opname,'_merge.txt',sep=''),ncol=ncol(mrgdat),append=T)

	## ATom GV+ 6.1 misses 7,925 because ObsPack has test flights that are not in merge product

} # end of function
