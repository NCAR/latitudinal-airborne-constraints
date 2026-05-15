# program to read in airborne trimmed 10-second merge data from observations or models, filter for strat, BL, and flight, further trim and output

### reading in published BL filt file?

filtac=function(oco2mipoutput=F,inclpfp=F,exclstrat=T,species='CO2_OP',useallrf=F,outdir='OBS',adj2X2019=F,tomfilenamein='RData.ATom.Mer.all_10X_trimmed',indir='OBS'){

	# define columns to keep (used by inclmed and inclpfp too)
	colsel=c('Year','DOY','UTC','PSXC','THETA','camp','flt','slc','CO2_AO2','CO2.X','CH4_NOAA','CH4_QCLS','N2O_QCLS','POT','UTC_Start','RF','G_LAT','G_LONG','G_ALT','P','A.no','O3_UCATS','O3_CL','H2O_UWV','H2O_DLH','H2Oppmv_vxl','N2O_P','N2O_PECD','N2O_UGC','N2O_UCATS','O3_ppb','O3_UO3','O3_UGC','VMR_VXL','H2O_NOAA','jday2016','UTC_Mean','YYYYMMDD','HDG','CO2_QCLS','CO2_NOAA')

	print(Sys.time());print('loading previously trimmed mergefile')
	print(tomfilenamein)
	if(!oco2mipoutput){ # observations
		load(paste(indir,'/',tomfilenamein,sep=''))
		if(grepl('CO2_MED',species)) atommerge=medmerge 
		if(grepl('CO2_PFP',species)) atommerge=pfpmerge 
		# > paste(names(atommerge),collapse=' ')
	} else { # OCO-2 MIP output
		atommerge=read.table(paste(indir,'/',tomfilenamein,sep=''),header=T)
		# year mon day hour min sec alt lat lon obspack_id camp flt prof pressure theta strat hdg co2
	}

	## Strat cutoff (uses N2O, O3, and H2O)
	if(exclstrat){

		print('filtering stratosphere')

		## strat flag is calculated in OBS/aircraft_obspack_merge.r and output to ATom_trimmed.txt and OBS/ATom_obspack_GLOBALVIEWplus_v6.1_merge.txt
		## strat flag is then merged with model output in V10MIP/model_aircraft_merge.r which reads in OBS/ATom_obspack_GLOBALVIEWplus_v6.1_merge.txt
		## for the obs, this program reads in the output from read_airborne_obs.r which has not yet had strat flag calculated, so recalculate here

		if(!oco2mipoutput){

			## based upon Jin et al., ACP, 2021 (https://doi.org/10.5194/acp-21-217-2021)
			stratcoh2o=50
			stratcoo3=150
			stratcon2o=319
			# for detrending n2o
			glbn2ofile='LOCALDATA/n2o_annmean_gl.txt' # from: https://gml.noaa.gov/webdata/ccgg/trends/n2o/n2o_annmean_gl.txt
			temp=readLines(glbn2ofile)
			hlines<-which(temp=='# year     mean      unc')
			glbn2o=read.table(glbn2ofile,skip=hlines,header=F,stringsAsFactors=F)
			colnames(glbn2o)=c('year','n2o','unc')

			atommerge$strat=rep(0,nrow(atommerge)) # 0 means trop
			h2oref=atommerge$H2O_DLH; h2oref[is.na(h2oref)]=atommerge$H2O_UWV[is.na(h2oref)]; h2oref[is.na(h2oref)]=0 # if H2O missing treat as if potentially strat
			n2oref=atommerge$N2O_QCLS; n2oref[is.na(n2oref)]=atommerge$N2O_PECD[is.na(n2oref)]; n2oref[is.na(n2oref)]=atommerge$N2O_UCATS[is.na(n2oref)]
			n2oref=n2oref-(approx(glbn2o$year+0.5,glbn2o$n2o,atommerge$Year+atommerge$DOY/365)$y-glbn2o$n2o[glbn2o$year==2009])
			n2oref[is.na(n2oref)]=400 # if N2O missing do not use for filter
			o3ref=atommerge$O3_CL; o3ref[is.na(o3ref)]=atommerge$O3_UCATS[is.na(o3ref)]; o3ref[is.na(o3ref)]=0 # if O3 missing do not use for filter
			atommerge$strat[h2oref<stratcoh2o&(o3ref>stratcoo3|n2oref<stratcon2o|(o3ref==0&n2oref==400))]=1 # if either o3 or n2o criteria are met, or if both are missing, consider strat
			atommerge$strat[h2oref==0&o3ref==0&n2oref==400&atommerge$G_ALT<8000]=0 # if all 3 missing assume < 8 km is trop
			tomstratfilt=atommerge$strat==0 # convert to old format (0/F = strat, 1/T = trop)

		} else {

			# already calculated in OBS/aircraft_obspack_merge.r and merged into model file in V10MIP/model_aircraft_merge.r
			tomstratfilt=atommerge$strat==0 # just reverse here (0/F = strat, 1/T = trop)
			
		}

		print(round(sum(!tomstratfilt&!is.na(tomstratfilt))/length(tomstratfilt),2))

	} else { # exclstrat=F

		tomstratfilt=rep(T,nrow(atommerge))

	}


	## Local / BL influence filter

	print(Sys.time());print('filtering BL over land')
	# filter out BL data over land (ANC, Fairbanks, CHC, Lauder)

	ints=read.csv('LOCALDATA/atom_filter_times.csv')
	# startdate            stopdate         reason
	# 1 2016-07-29 14:32:18 2016-07-29 14:34:48        takeoff
	intstartdt=strptime(ints$startdate,format='%Y-%m-%d %H:%M:%S',tz='UTC')
	intstopdt=strptime(ints$stopdate,format='%Y-%m-%d %H:%M:%S',tz='UTC')
	tomblfilt=rep(T,nrow(atommerge))

	if(!oco2mipoutput){
		tomdt=strptime(paste(atommerge$Year,'-',atommerge$DOY,sep=''),format='%Y-%j',tz='UTC')+atommerge$UTC # DOY is day of take-off and UTC is seconds since midnight before take-off
	} else {
		tomdt=ISOdatetime(atommerge$year,atommerge$mon,atommerge$day,atommerge$hour,atommerge$min,atommerge$sec,tz='UTC')
	}

	for(i in c(1:nrow(ints))){
		tomblfilt[difftime(tomdt,intstartdt[i])>=0&difftime(tomdt,intstopdt[i])<=0]=F
	}

	print(paste('ATom BL filt removing',sum(!tomblfilt),'out of',length(tomblfilt),',',round(sum(!tomblfilt)/length(tomblfilt)*100,1),'%'))


	# flights to include
	print(Sys.time());print('selecting flights')
	if(useallrf){

		tomfltsel=atommerge$flt>0

	} else { # restrict to flights used in xsect only

		tomfltsel=rep(F,nrow(atommerge))
		if(oco2mipoutput) atommerge$A.no=atommerge$camp
		anum=1
		fnum=c(1:9) ## excluding 10 and 11 across N. America
		tomfltsel[atommerge$A.no==anum&is.element(atommerge$flt,fnum)]=T
		anum=2
		fnum=c(1:11) ## excluding ANC- ferry
		tomfltsel[atommerge$A.no==anum&is.element(atommerge$flt,fnum)]=T
		anum=3
		fnum=c(1:13) ## excluding ANC- ferry
		tomfltsel[atommerge$A.no==anum&is.element(atommerge$flt,fnum)]=T
		anum=4
		fnum=c(1:13) ## excluding ANC- ferry
		tomfltsel[atommerge$A.no==anum&is.element(atommerge$flt,fnum)]=T

	}

	atommerge=atommerge[tomfltsel&tomblfilt&tomstratfilt,]
	print(dim(atommerge))

	# make common campaign variable
	if(!oco2mipoutput){
		atommerge$camp=atommerge$A.no+20
	} else {
		atommerge$camp=atommerge$camp+20
	}

	# assign slice by heading (0 = Southbound, 1 = Northbound)
	if(oco2mipoutput) atommerge$HDG=atommerge$hdg
	atommerge$slc=(atommerge$HDG<90|atommerge$HDG>=270)*1

	# trim and then merge
	if(!oco2mipoutput){
		tomtrim=atommerge[,c('Year','DOY','UTC','G_LAT','P','camp','flt','slc',species)]; colnames(tomtrim)=c('Year','DOY','UTC','GGLAT','PSXC','camp','flt','slc',species)
	} else {
		tomtrim=atommerge
		tomtrim$GGLAT=tomtrim$lat
		tomtrim$PSXC=tomtrim$pressure
		tomtrim$CO2_OP=tomtrim$co2
	}

	mergefile=tomtrim

	save('mergefile',file=paste(outdir,'/filtcombmergeout_',species,sep=''))

	return(mergefile)

} # end of filtac() function
