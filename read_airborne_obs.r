# program to read in airborne 10-second merge data from models, trim for selected variables and output


# note on scales: ATom observations and OCO-2 V10MIP output are on WMO X2007, so adj2X2019 set to F here and adjMLO2X2007 set to T in harm_fit.r (for MLO). 
# Medusa CO2 data also require adjustment from SIO CO2 to WMO X2007 scale (after which adj2X2019 determines if they are also adjusted to X2019)

readac=function(adj2X2019=F,inclmed=T,inclpfp=T,outdir='OBS',tomfilenameout='RData.ATom.Mer.all_10X_trimmed',medfilenameout='RData.ATom.Mer.all_MED_trimmed',pfpfilenameout='RData.ATom.Mer.all_PFP_trimmed'){

	# define columns to keep (used by inclmed and inclpfp too)
	colsel=c('A.no','camp','slice','RF','flt','Year','Month','Day','Hour','Min','Sec','DOY','jday2016','YYYYMMDD','UTC_Start','UTC_Mean','UTC','G_LAT','G_LONG','G_ALT','HDG','P','PSXC','THETA','POT','CO2_NOAA','CO2_QCLS','CO2_AO2','CO2.X','CH4_NOAA','CH4_QCLS','N2O_QCLS','N2O_P','N2O_PECD','N2O_UGC','N2O_UCATS','H2O_UWV','H2O_DLH','H2Oppmv_vxl','VMR_VXL','H2O_NOAA','O3_UCATS','O3_CL','O3_ppb','O3_UO3','O3_UGC')

	print(Sys.time());print('reading in mergefiles')

	## read in ATom files
	library('ncdf4')
	atommergedir='LOCALDATA/ATOM' # these are version 2.0 (21-08-26)
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
	atommerge$Year=as.numeric(substr(atommerge$Flight_Date,1,4)) # Flight_Date is day of takeoff
	atommerge$Month=as.numeric(substr(atommerge$Flight_Date,5,6)) # Flight_Date is day of takeoff
	atommerge$Day=as.numeric(substr(atommerge$Flight_Date,7,8)) # Flight_Date is day of takeoff
	atommerge$UTC_Mean=atommerge$UTC_Start+5 # UTC_Start is seconds since midnight prior to takeoff (values greater than 86400 when flight goes over the following midnight)
	atomdt=as.POSIXlt(ISOdatetime(atommerge$Year,atommerge$Month,atommerge$Day,0,0,0,tz='UTC')+atommerge$UTC_Mean,tz='UTC') # use day of takeoff here since adding seconds since midnight on day of takeoff
	# reset to calculate actual UTC Day for portion of flights over midnight
	atommerge$Year=as.POSIXlt(atomdt)$year+1900  
	atommerge$Month=as.POSIXlt(atomdt)$mon+1
	atommerge$Day=as.POSIXlt(atomdt)$mday # actual day (steps over midnight)
	atommerge$Hour=as.POSIXlt(atomdt)$hour
	atommerge$Min=as.POSIXlt(atomdt)$min
	atommerge$Sec=as.POSIXlt(atomdt)$sec
	atommerge$DOY=as.POSIXlt(atomdt)$yday+1 # Jan. 1 = 1
	atommerge$DOY[atommerge$UTC_Mean>86400]=atommerge$DOY[atommerge$UTC_Mean>86400]-1 ## DOY is day of year of takeoff
	atommerge$jday2016=atommerge$time/86400
	atommerge$YYYYMMDD=strftime(as.POSIXlt(atomdt),format='%Y%m%d',tz='UTC')

	atommerge=data.frame(atommerge)

	atommerge$THETA=atommerge$POT
	atommerge$UTC=trunc(atommerge$UTC_Start)
	atommerge$flt=atommerge$RF
	atommerge$PSXC=atommerge$P

	if(adj2X2019){ ## convert aircraft obs from X2007 to X2019 scale
		# from NOAA webpage: X2019 = 1.00079 * X2007 - 0.142
		atommerge$CO2_NOAA=atommerge$CO2_NOAA*1.00079-0.142
		atommerge$CO2_QCLS=atommerge$CO2_QCLS*1.00079-0.142
		atommerge$CO2_AO2=atommerge$CO2_AO2*1.00079-0.142
	}

	## create ATom CO2 difference variables
	atommerge$CO2_AO2_m_CO2_NOAA=atommerge$CO2_AO2-atommerge$CO2_NOAA
	atommerge$CO2_QCLS_m_CO2_NOAA=atommerge$CO2_QCLS-atommerge$CO2_NOAA

	## make CO2_OP
	atommerge$CO2_OP=atommerge$CO2_NOAA # matching obspack
	atommerge$CO2_OP[atommerge$A.no==1&(atommerge$RF==1|atommerge$RF==2)]=atommerge$CO2_QCLS[atommerge$A.no==1&(atommerge$RF==1|atommerge$RF==2)]

	colsel=c(colsel,'CO2_AO2_m_CO2_NOAA','CO2_QCLS_m_CO2_NOAA','CO2_OP')

	atommerge=atommerge[,match(colsel,names(atommerge),nomatch=F)] # reorder columns
	save('atommerge',file=paste(outdir,'/',tomfilenameout,sep=''))

	if(inclmed){

		# add in MEDUSA data
		print(Sys.time());print('reading in Medusa')

		colsel=c(colsel,'CO2_MED','CO2_MED_m_CO2_NOAA')
		# MED files also include in situ tracers, these get added inline to 10-sec variables - leave to allow for difference calculations

		if(T){

			medmerge=NULL
			for(anum in c(1:4)){

				print(paste('A',anum,sep=''))
				campaign=paste('ATOM',anum,sep='')
				campdata=read.table(paste('LOCALDATA/MEDUSA_MERGE_',campaign,'_200918.tbl',sep=''),header=T)
				campdata$A.no=rep(anum,nrow(campdata))
				medmerge=rbind(medmerge,campdata)

			} # loop on anum

			## to adjust, use 2015 NCAR primary comparison (most recent)
			# from /h/eol/stephens/NCAR/CALFAC/PRIMARIES/1912bkup/gmd-sio.r.Rout
			# [1] "2015 fit:"
			# [1]  1.082731e+00 -1.019421e-02  1.795693e-05
			sio2gmd2=-1.795693e-05; sio2gmd1=1.019421e-02; sio2gmd0=-1.082731e+00
			# print(400^2*sio2gmd2+400*sio2gmd1+sio2gmd0)
			# at 400 ppm, Scripps O2 lab low by 0.12 ppm
			print('Adjusting ATom CO2_MED to WMO (X2007) Scale, based on 2015 SIO/GMD analysis of NCAR primaries')
			medmerge$CO2_MED=medmerge$CO2_MED^2*sio2gmd2+medmerge$CO2_MED*sio2gmd1+sio2gmd0+medmerge$CO2_MED
			print(paste('offset at 400 ppm =',round(400^2*sio2gmd2+400*sio2gmd1+sio2gmd0,3)))

			if(adj2X2019){ ## convert aircraft obs to X2019 scale
				# from NOAA webpage: X2019 = 1.00079 * X2007 - 0.142
				medmerge$CO2_MED=medmerge$CO2_MED*1.00079-0.142
				medmerge$CO2_NOAA=medmerge$CO2_NOAA*1.00079-0.142
				medmerge$CO2_QCLS=medmerge$CO2_QCLS*1.00079-0.142
				medmerge$CO2_AO2=medmerge$CO2_AO2*1.00079-0.142
			}
			# create ATom CO2 difference variable
			medmerge$CO2_NOAA[medmerge$CO2_NOAA>1000]=NA #### two values 1962.9380 26723.8600 for some reason (rows 159:160 of ATOM1 file)
			medmerge$CO2_MED_m_CO2_NOAA=medmerge$CO2_MED-medmerge$CO2_NOAA
			medmerge$G_ALT=medmerge$GGALT
			medmerge$G_LAT=medmerge$GGLAT
			medmerge$G_LON=medmerge$GGLON
			## DOY in LOCALDATA/MEDUSA_MERGE_',campaign,'_200918.tbl is messed up, have to regenerate
			medmerge$DOY=trunc(medmerge$jday2016)+1 # jday2016 is day of year since Dec. 31 2015 (Jan. 1, 2016 = 1) and fraction of year
			medmerge$DOY[medmerge$A.no==2|medmerge$A.no==3]=trunc(medmerge$jday2016[medmerge$A.no==2|medmerge$A.no==3])+1-366
			medmerge$DOY[medmerge$A.no==4]=trunc(medmerge$jday2016[medmerge$A.no==4])+1-366-365
			medmerge$DOY[medmerge$UTC_Mean>86400]=medmerge$DOY[medmerge$UTC_Mean>86400]-1 ## do not want DOY to jump a whole day over midnight
			#print(paste('ATom-',anum,' CO2_AO2-CO2_MED:',round(mean(medmerge$CO2_AO2-medmerge$CO2_MED,na.rm=T),3),sep=''))

			medmerge=data.frame(medmerge,stringsAsFactors=F)

		} else { # code to use public version of Medusa merge files. Those files have merge errors, so using internal (correct) version here.

			medmerge=NULL
			medmergefiles=c('MER-MED_DC8_ATom-1.nc','MER-MED_DC8_ATom-2.nc','MER-MED_DC8_ATom-3.nc','MER-MED_DC8_ATom-4.nc') 
			#medvar=c(atomvar,'MED/CO2_MED')
			medvar=c('UTC_Start','Flight_Date','time','DLH-H2O/H2O_ppmv','UCATS-H2O/H2O_UWV','QCLS-CH4-CO-N2O/N2O_QCLS','NOyO3-O3/O3_CL','UCATS-O3/O3_UO3','MMS/G_ALT','RF','prof.no','MMS/P','MMS/POT','MMS/G_LAT','MMS/G_LONG','NOAA-Picarro/CO2_NOAA','NOAA-Picarro.wt/CO2_NOAA.wt','QCLS-CO2/CO2_QCLS','AO2/CO2_AO2','NOAA-Picarro/CH4_NOAA','QCLS-CH4-CO-N2O/CH4_QCLS','MEDUSA/CO2_MED') # no GCECD, UCATS-GC, MMS/HDG, or CO2.X. DLH-H2O/H2O_ppmv and UCATS-O3_UO3 different names
			### currently crashing filter_airborne.r which looks for N2O_PECD
			for(i in c(1:4)){
				mednc=nc_open(paste(atommergedir,'/',medmergefiles[i],sep=''))
				count=length(ncvar_get(mednc,'UTC_Start'))
				campdata=NULL
				for(var in medvar){
					if(i==1&var=='UCATS-H2O/H2O_UWV'){ # no UCATS H2O on ATom-1
						campdata=cbind(campdata,rep(NA,count))
					} else {
						campdata=cbind(campdata,ncvar_get(mednc,var))
					}
				}
				nc_close(mednc)
				campdata=cbind(campdata,rep(i,count)) # A.no
				medmerge=rbind(medmerge,campdata)
			}
			medmerge=data.frame(medmerge,stringsAsFactors=F)
			names(medmerge)=c(gsub('.*/','',medvar),'A.no')

			## Adjust MED_CO2 from SIO to WMO X2007 scale
			## to adjust, use 2015 NCAR primary comparison (most recent)
			# from /h/eol/stephens/NCAR/CALFAC/PRIMARIES/1912bkup/gmd-sio.r.Rout
			# [1] "2015 fit:"
			# [1]  1.082731e+00 -1.019421e-02  1.795693e-05
			sio2gmd2=-1.795693e-05; sio2gmd1=1.019421e-02; sio2gmd0=-1.082731e+00
			# print(400^2*sio2gmd2+400*sio2gmd1+sio2gmd0)
			# at 400 ppm, Scripps O2 lab low by 0.12 ppm
			print('Adjusting ATom CO2_MED to WMO (X2007) Scale, based on 2015 SIO/GMD analysis of NCAR primaries')
			medmerge$CO2_MED=medmerge$CO2_MED^2*sio2gmd2+medmerge$CO2_MED*sio2gmd1+sio2gmd0+medmerge$CO2_MED
			print(paste('offset at 400 ppm =',round(400^2*sio2gmd2+400*sio2gmd1+sio2gmd0,3)))

			if(adj2X2019){ ## convert aircraft obs to X2019 scale
				# from NOAA webpage: X2019 = 1.00079 * X2007 - 0.142
				medmerge$CO2_MED=medmerge$CO2_MED*1.00079-0.142
				medmerge$CO2_NOAA=medmerge$CO2_NOAA*1.00079-0.142
				medmerge$CO2_QCLS=medmerge$CO2_QCLS*1.00079-0.142
				medmerge$CO2_AO2=medmerge$CO2_AO2*1.00079-0.142
			}

			# add time variables
			medmerge$Year=as.numeric(substr(medmerge$Flight_Date,1,4)) # Flight_Date is day of takeoff
			medmerge$Month=as.numeric(substr(medmerge$Flight_Date,5,6)) # Flight_Date is day of takeoff
			medmerge$Day=as.numeric(substr(medmerge$Flight_Date,7,8)) # Flight_Date is day of takeoff
			medmerge$UTC_Mean=medmerge$UTC_Start+5 # UTC_Start is seconds since midnight prior to takeoff (values greater than 86400 when flight goes over the following midnight)
			meddt=as.POSIXlt(ISOdatetime(medmerge$Year,medmerge$Month,medmerge$Day,0,0,0,tz='UTC')+medmerge$UTC_Mean,tz='UTC') # use day of takeoff here since adding seconds since midnight on day of takeoff
			# reset to calculate actual UTC Day for portion of flights over midnight
			medmerge$Year=as.POSIXlt(meddt)$year+1900
			medmerge$Month=as.POSIXlt(meddt)$mon+1
			medmerge$Day=as.POSIXlt(meddt)$mday # actual day (steps over midnight)
			medmerge$Hour=as.POSIXlt(meddt)$hour
			medmerge$Min=as.POSIXlt(meddt)$min
			medmerge$Sec=as.POSIXlt(meddt)$sec
			medmerge$DOY=as.POSIXlt(meddt)$yday+1 # Jan. 1 = 1
			medmerge$DOY[medmerge$UTC_Mean>86400]=medmerge$DOY[medmerge$UTC_Mean>86400]-1 ## DOY is day of year of takeoff
			medmerge$jday2016=medmerge$time/86400
			medmerge$YYYYMMDD=strftime(as.POSIXlt(meddt),format='%Y%m%d',tz='UTC')

			medmerge$THETA=medmerge$POT
			medmerge$UTC=trunc(medmerge$UTC_Start)
			medmerge$flt=medmerge$RF
			medmerge$PSXC=medmerge$P
			medmerge$THETA=medmerge$POT

			if(adj2X2019){ ## convert aircraft obs from X2007 to X2019 scale
				# from NOAA webpage: X2019 = 1.00079 * X2007 - 0.142
				medmerge$CO2_NOAA=medmerge$CO2_NOAA*1.00079-0.142
				medmerge$CO2_QCLS=medmerge$CO2_QCLS*1.00079-0.142
				medmerge$CO2_AO2=medmerge$CO2_AO2*1.00079-0.142
			}

			## create ATom CO2 difference variables
			medmerge$CO2_MED_m_CO2_NOAA=medmerge$CO2_MED-medmerge$CO2_NOAA

		} # if(F)

		medmerge=medmerge[,match(colsel,names(medmerge),nomatch=F)] # reorder columns
		save('medmerge',file=paste(outdir,'/',medfilenameout,sep=''))

		print(c(sum(!is.na(medmerge$CO2_MED)),sum(!is.na(medmerge$CO2_NOAA)),sum(!is.na(medmerge$CO2_MED_m_CO2_NOAA))))

	} # if(inclmed)


	if(inclpfp){
		# add in PFP data
		print(Sys.time());print('reading in PFP')

		colsel=c(colsel,'CO2_PFP','CO2_PFP_m_CO2_NOAA')
		# PFP files also include in situ tracers, these get added inline to 10-sec variables - leave to allow for difference calculations

		if(F){

			pfp=read.csv('LOCALDATA/Mor.PFP.all.at1234.2021-06-13.tbl')
			### read public file or move locally

			pfp$DOY=trunc(pfp$jday2016)+1 # jday2016 is day of year since Dec. 31 2015 (Jan. 1, 2016 = 1) and fraction of year
			pfp$DOY[pfp$A.no==2|pfp$A.no==3]=trunc(pfp$jday2016[pfp$A.no==2|pfp$A.no==3])+1-366
			pfp$DOY[pfp$A.no==4]=trunc(pfp$jday2016[pfp$A.no==4])+1-366-365
			pfp$DOY[pfp$UTC_Mean>86400]=pfp$DOY[pfp$UTC_Mean>86400]-1 ## do not want DOY to jump a whole day over midnight
			pfp$Year=as.numeric(substr(pfp$YYYYMMDD,1,4))

			pfp$THETA=pfp$POT
			pfp$UTC=trunc(pfp$UTC_Start)
			pfp$flt=pfp$RF
			pfp$PSXC=pfp$P
				
			if(adj2X2019){ ## convert aircraft obs to X2019 scale
				# from NOAA webpage: X2019 = 1.00079 * X2007 - 0.142
				pfp$CO2_PFP=pfp$CO2_PFP*1.00079-0.142
				pfp$CO2_NOAA=pfp$CO2_NOAA*1.00079-0.142
				pfp$CO2_QCLS=pfp$CO2_QCLS*1.00079-0.142
				pfp$CO2_AO2=pfp$CO2_AO2*1.00079-0.142
			}

			# create ATom CO2 difference variable
			pfp$CO2_PFP_m_CO2_NOAA=pfp$CO2_PFP-pfp$CO2_NOAA

			# trim to speed up:
			pfp=pfp[,is.element(colnames(pfp),colsel)]

		} else { # code to use public version of PFP merge files

			pfpmerge=NULL
			pfpmergefiles=c('MER-PFP_DC8_ATom-1.nc','MER-PFP_DC8_ATom-2.nc','MER-PFP_DC8_ATom-3.nc','MER-PFP_DC8_ATom-4.nc') 
			pfpvar=c(atomvar,'PFP/CO2_PFP')
#			pfpvar=c('UTC_Start','Flight_Date','time','DLH-H2O/H2O_ppmv','UCATS-H2O/H2O_UWV','QCLS-CH4-CO-N2O/N2O_QCLS','NOyO3-O3/O3_CL','UCATS-O3/O3_UO3','MMS/G_ALT','RF','prof.no','MMS/P','MMS/POT','MMS/G_LAT','MMS/G_LONG','NOAA-Picarro/CO2_NOAA','NOAA-Picarro.wt/CO2_NOAA.wt','QCLS-CO2/CO2_QCLS','AO2/CO2_AO2','NOAA-Picarro/CH4_NOAA','QCLS-CH4-CO-N2O/CH4_QCLS','PFP/CO2_PFP') # no GCECD, UCATS-GC, MMS/HDG, or CO2.X. DLH-H2O/H2O_ppmv and UCATS-O3_UO3 different names
			for(i in c(1:4)){
				pfpnc=nc_open(paste(atommergedir,'/',pfpmergefiles[i],sep=''))
				count=length(ncvar_get(pfpnc,'UTC_Start'))
				campdata=NULL
				for(var in pfpvar){
					if(i==1&var=='UCATS-H2O/H2O_UWV'){ # no UCATS H2O on ATom-1
						campdata=cbind(campdata,rep(NA,count))
					} else {
						campdata=cbind(campdata,ncvar_get(pfpnc,var))
					}
				}
				nc_close(pfpnc)
				campdata=cbind(campdata,rep(i,count)) # A.no
				pfpmerge=rbind(pfpmerge,campdata)
			}
			pfpmerge=data.frame(pfpmerge,stringsAsFactors=F)
			names(pfpmerge)=c(gsub('.*/','',pfpvar),'A.no')

			if(adj2X2019){ ## convert aircraft obs to X2019 scale
				# from NOAA webpage: X2019 = 1.00079 * X2007 - 0.142
				pfpmerge$CO2_PFP=pfpmerge$CO2_PFP*1.00079-0.142
				pfpmerge$CO2_NOAA=pfpmerge$CO2_NOAA*1.00079-0.142
				pfpmerge$CO2_QCLS=pfpmerge$CO2_QCLS*1.00079-0.142
				pfpmerge$CO2_AO2=pfpmerge$CO2_AO2*1.00079-0.142
			}

			# add time variables
			pfpmerge$Year=as.numeric(substr(pfpmerge$Flight_Date,1,4)) # Flight_Date is day of takeoff
			pfpmerge$Month=as.numeric(substr(pfpmerge$Flight_Date,5,6)) # Flight_Date is day of takeoff
			pfpmerge$Day=as.numeric(substr(pfpmerge$Flight_Date,7,8)) # Flight_Date is day of takeoff
			pfpmerge$UTC_Mean=pfpmerge$UTC_Start+5 # UTC_Start is seconds since midnight prior to takeoff (values greater than 86400 when flight goes over the following midnight)
			pfpdt=as.POSIXlt(ISOdatetime(pfpmerge$Year,pfpmerge$Month,pfpmerge$Day,0,0,0,tz='UTC')+pfpmerge$UTC_Mean,tz='UTC') # use day of takeoff here since adding seconds since midnight on day of takeoff
			# reset to calculate actual UTC Day for portion of flights over midnight
			pfpmerge$Year=as.POSIXlt(pfpdt)$year+1900
			pfpmerge$Month=as.POSIXlt(pfpdt)$mon+1
			pfpmerge$Day=as.POSIXlt(pfpdt)$mday # actual day (steps over midnight)
			pfpmerge$Hour=as.POSIXlt(pfpdt)$hour
			pfpmerge$Min=as.POSIXlt(pfpdt)$min
			pfpmerge$Sec=as.POSIXlt(pfpdt)$sec
			pfpmerge$DOY=as.POSIXlt(pfpdt)$yday+1 # Jan. 1 = 1
			pfpmerge$DOY[pfpmerge$UTC_Mean>86400]=pfpmerge$DOY[pfpmerge$UTC_Mean>86400]-1 ## DOY is day of year of takeoff
			pfpmerge$jday2016=pfpmerge$time/86400
			pfpmerge$YYYYMMDD=strftime(as.POSIXlt(pfpdt),format='%Y%m%d',tz='UTC')

			pfpmerge$THETA=pfpmerge$POT
			pfpmerge$UTC=trunc(pfpmerge$UTC_Start)
			pfpmerge$flt=pfpmerge$RF
			pfpmerge$PSXC=pfpmerge$P

			## create ATom CO2 difference variables
			pfpmerge$CO2_PFP_m_CO2_NOAA=pfpmerge$CO2_PFP-pfpmerge$CO2_NOAA

		} # if(T)

		pfpmerge=pfpmerge[,match(colsel,names(pfpmerge),nomatch=F)] # reorder columns
		save('pfpmerge',file=paste(outdir,'/',pfpfilenameout,sep=''))

	} # if(inclpfp)

} # end of readac() function
