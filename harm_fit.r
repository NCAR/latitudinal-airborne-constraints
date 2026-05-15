# program to read in trimmed and filtered airborne 10-second observations or models, bin by pressure and latitude, calculate harmonic fits, and output results
## also option to detrend

# note on scales: ATom V2 merge and OCO-2 V10MIP output are on WMO X2007, so adj2X2019 set to F in read_airborne.r and adjMLO2X2007 set to T here
# Medusa CO2 data also require adjustment from SIO CO2 to WMO X2007 scale (after which adj2X2019 determines if they are also adjusted to X2019)

source('filter_airborne.r')

harmf=function(nharm=1,oco2mipoutput=F,prsspan=100,prsrange=c(100,1000),latrange=c(-70,90),latspan=10,calatco=20,species='APO_AO2',surfco=900,coltop=300,colwtco=0.7,detrend=T,dtrsta='mlo',flag='',aggby='campflt',ylim1='byalt',units='per meg',binxaxcex=1.0,binyaxcex=1.0,binprscex=1.0,binlatcex=1.0,outdir='OBS',reqncamp=T,indir='OBS',tomfilenamein='RData.ATom.Mer.all_10X_trimmed'){

	# reqncamp = whether to require the number of campaigns in a bin to be > nharm*2

	print(Sys.time());print('filtering and combining previously trimmed mergefile')
	mergefile=filtac(species=species,indir=indir,outdir=outdir,tomfilenamein=tomfilenamein,oco2mipoutput=oco2mipoutput)

	if(!oco2mipoutput){
		mergefiledt=as.POSIXlt(strptime(paste(mergefile$Year,' ',mergefile$DOY,' 00:00:00',sep='',tz='UTC'),format='%Y %j %H:%M:%S')+mergefile$UTC,tz='UTC')
	} else {
		mergefiledt=as.POSIXlt(ISOdatetime(mergefile$year,mergefile$mon,mergefile$day,mergefile$hour,mergefile$min,mergefile$sec,tz='UTC'))
	}

	if(detrend){

		print(Sys.time());print('detrending')

		# read in output from noaa_stl.r for detrending
		statrd=read.table(paste(dtrsta,'_co2_stl_results.txt',sep=''),header=T,stringsAsFactors=F)

		stadt=strptime(paste('16',statrd$date,sep=''),format='%d%b%Y',tz='UTC')
		statrd$year=stadt$year+1900
		statrd$mon=stadt$mon+1
		statrd$ydays=rep(365,length(statrd$year))
		statrd$ydays[statrd$year%%4==0]=366
		statrd$yrfrac=statrd$year+stadt$yday/statrd$ydays
		statrd$ti=statrd$trend+statrd$iav
		trendcorr=approx(as.POSIXct(stadt),statrd$ti,as.POSIXct(mergefiledt),rule=1)$y
		if(adjMLO2X2007){ ## convert MLO back to X2007 scale
			# from NOAA webpage: X2019 = 1.00079 * X2007 - 0.142
			# thus, X2007 = (X2019 + 0.142) / 1.00079
			trendcorr=(trendcorr+0.142)/1.00079 ## minus ~ 0.2 at 400 ppm
		}
		mergefile[,species]=mergefile[,species]-trendcorr

	}

	print(Sys.time());print('aggregating')
	if(length(latspan)==1){
		latedges=seq(latrange[1],latrange[2],latspan)
	} else {
		latedges=latspan
	}
	latvals=filter(latedges,c(0.5,0.5));latvals=latvals[!is.na(latvals)]
	latmin=latedges[1:(length(latedges)-1)]
	latmax=latedges[2:length(latedges)]
	nlat=length(latvals)
	prsedges=seq(prsrange[1],prsrange[2],prsspan)
	prsvals=filter(prsedges,c(0.5,0.5));prsvals=prsvals[!is.na(prsvals)]
	prsmin=prsedges[1:(length(prsedges)-1)]
	prsmax=prsedges[2:length(prsedges)]
	nprs=length(prsvals)
	latbin=cut(mergefile$GGLAT,latedges,labels=F)
	prsbin=cut(mergefile$PSXC,prsedges,labels=F)

	if(aggby=='day'){
		binval=aggregate(mergefile[,species],by=list(lat=latbin,prs=prsbin,yday=mergefiledt$yday),mean,na.rm=T)
		bincount=aggregate(!is.na(mergefile[,species]),by=list(lat=latbin,prs=prsbin,yday=mergefiledt$yday),sum,na.rm=T)
	} else if(aggby=='campflt'){
		binval=aggregate(mergefile[,species],by=list(lat=latbin,prs=prsbin,camp=mergefile$camp,flt=mergefile$flt),mean,na.rm=T)
		bincount=aggregate(!is.na(mergefile[,species]),by=list(lat=latbin,prs=prsbin,camp=mergefile$camp,flt=mergefile$flt),sum,na.rm=T)
		binval$yday=aggregate(as.numeric(mergefiledt$yday),by=list(lat=latbin,prs=prsbin,camp=mergefile$camp,flt=mergefile$flt),mean,na.rm=T)$x
	} else if(aggby=='campslc'){
		binval=aggregate(mergefile[,species],by=list(lat=latbin,prs=prsbin,camp=mergefile$camp,slc=mergefile$slc),mean,na.rm=T)
		bincount=aggregate(!is.na(mergefile[,species]),by=list(lat=latbin,prs=prsbin,camp=mergefile$camp,slc=mergefile$slc),sum,na.rm=T)
		binval$yday=aggregate(as.numeric(mergefiledt$yday),by=list(lat=latbin,prs=prsbin,camp=mergefile$camp,slc=mergefile$slc),mean,na.rm=T)$x
	}

	write(colnames(binval),paste(outdir,'/Bin_Values_',species,flag,'.txt',sep=''),ncol=ncol(binval))
	write(rbind(latvals[binval$lat],prsvals[binval$prs],t(binval[,3:ncol(binval)])),paste(outdir,'/Bin_Values_',species,flag,'.txt',sep=''),append=T,ncol=ncol(binval))

	print(Sys.time());print('fitting bin harmonics')
	# generage harmonics
	yrfrac=(binval$yday+0.5)/365*2*pi
	harm=cbind(cos(yrfrac),sin(yrfrac),cos(2*yrfrac),sin(2*yrfrac),cos(3*yrfrac),sin(3*yrfrac),cos(4*yrfrac),sin(4*yrfrac))[,1:(nharm*2)]
	xdays=seq(0.5,364.5)/365*2*pi # noon every day of non-leap year
	hdays=cbind(rep(1,length(xdays)),cos(xdays),sin(xdays),cos(2*xdays),sin(2*xdays),cos(3*xdays),sin(3*xdays),cos(4*xdays),sin(4*xdays))[,1:(nharm*2+1)]

	pngw=1800;pngh=1200;psize=24
	png(paste(outdir,'/Bin_Fits_',species,flag,'.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	par(mfrow=c(nprs,nlat))
	par(mar=c(0.1,0.1,0.1,0.1))
	par(oma=c(4,4,4,4))
	par(mgp=c(2.5,0.75,0))

	seasamp=matrix(NA,nlat,nprs)
	annmean=matrix(NA,nlat,nprs)
	seasphs=matrix(NA,nlat,nprs)
	count=matrix(NA,nlat,nprs)
	alldaily=array(NA,c(nlat,nprs,365))

	write(paste('lat prs annmn',paste(c('cos','sin','cos2','sin2','cos3','sin3','cos4','sin4')[1:(nharm*2)],collapse=' ')),paste(outdir,'/Fit_Coefficients_',species,flag,'.txt',sep=''))

	# do once to calculate y-range
	for(prs in c(1:nprs)){
		for(lat in c(1:nlat)){
			numpts=sum(!is.na(binval$x[binval$lat==lat&binval$prs==prs]))
			numcamp=length(unique(binval$camp[binval$lat==lat&binval$prs==prs&!is.na(binval$x)]))
			if(numpts>nharm*2&(!reqncamp|numcamp>nharm*2)){
				# fit harmonics
				coefs=lm(binval$x[binval$lat==lat&binval$prs==prs] ~ harm[binval$lat==lat&binval$prs==prs,])$coef
				# regenerate harmonics
				daily=drop(coefs%*%t(hdays))
				alldaily[lat,prs,]=daily
				seasamp[lat,prs]=diff(range(daily))
				annmean[lat,prs]=coefs[1]
				seasphs[lat,prs]=seq(1,365)[daily==min(daily)][1] # in case of ties, pick first one
				count[lat,prs]=sum(bincount$x[bincount$lat==lat&bincount$prs==prs])
				write(c(lat,prs,coefs),paste(outdir,'/Fit_Coefficients_',species,flag,'.txt',sep=''),ncol=length(coefs)+2,append=T)
			} else { 
				if(grepl('MED_m_',species)|grepl('PFP_m_',species)){ # for flask-insitu differences, just fit a constant
					binmean=mean(binval$x[binval$lat==lat&binval$prs==prs],na.rm=T)
					daily=rep(binmean,365)
					alldaily[lat,prs,]=daily
					seasamp[lat,prs]=0
					annmean[lat,prs]=binmean
					seasphs[lat,prs]=1
					count[lat,prs]=sum(bincount$x[bincount$lat==lat&bincount$prs==prs])
					write(c(lat,prs,rep(0,nharm*2),binmean),paste(outdir,'/Fit_Coefficients_',species,flag,'.txt',sep=''),ncol=nharm*2+3,append=T)
				} else {
					write(c(lat,prs,rep(NA,nharm*2+1)),paste(outdir,'/Fit_Coefficients_',species,flag,'.txt',sep=''),ncol=nharm*2+3,append=T)
				}
			}
		} # loop on lat
	} # loop on prs
	if(ylim1=='fixed') ylm=range(alldaily,na.rm=T)

	# redo for plotting
	for(prs in c(1:nprs)){
		if(ylim1=='byalt') ylm=range(c(alldaily[,prs,],binval$x[binval$prs==prs]),na.rm=T)
		for(lat in c(1:nlat)){
			numpts=sum(!is.na(binval$x[binval$lat==lat&binval$prs==prs]))
			numcamp=length(unique(binval$camp[binval$lat==lat&binval$prs==prs&!is.na(binval$x)]))
			if(numpts>nharm*2&(!reqncamp|numcamp>nharm*2)){
				if(ylim1=='bypanel') ylm=range(binval$x[binval$lat==lat&binval$prs==prs],na.rm=T)
				# fit harmonics
				coefs=lm(binval$x[binval$lat==lat&binval$prs==prs] ~ harm[binval$lat==lat&binval$prs==prs,])$coef
				# regenerate harmonics
				daily=drop(coefs%*%t(hdays))
				plot((binval$yday[binval$lat==lat&binval$prs==prs]+0.5)/365,binval$x[binval$lat==lat&binval$prs==prs],main='',xlab='',ylab='',axes=F,cex.main=0.75,cex.axis=0.25,cex.lab=1.00,ylim=ylm,xlim=c(0,1),type='n')
				lines(seq(0.5,364.5)/365,daily,col='black')
				points((binval$yday[binval$lat==lat&binval$prs==prs&binval$camp>20]+0.5)/365,binval$x[binval$lat==lat&binval$prs==prs&binval$camp>20],col='dark blue') # ATom
				# binval$camp values for ATom = 21 22 23 24
				box()
				if(lat==1&ylim1!='bypanel') axis(2,cex.axis=binyaxcex)
				if(prs==nprs) axis(1,cex.axis=binxaxcex)
			} else {
				if(grepl('MED_m_',species)|grepl('PFP_m_',species)){ # for flask-insitu differences, just fit a constant
					binmean=mean(binval$x[binval$lat==lat&binval$prs==prs],na.rm=T)
					daily=rep(binmean,365)
					plot((binval$yday[binval$lat==lat&binval$prs==prs]+0.5)/365,binval$x[binval$lat==lat&binval$prs==prs],main='',xlab='',ylab='',axes=F,cex.main=0.75,cex.axis=0.25,cex.lab=1.00,ylim=ylm,xlim=c(0,1),type='n')
					lines(seq(0.5,364.5)/365,daily,col='black')
					points((binval$yday[binval$lat==lat&binval$prs==prs&binval$camp>20]+0.5)/365,binval$x[binval$lat==lat&binval$prs==prs&binval$camp>20],col='dark blue') # ATom
					box()
					if(lat==1&ylim1!='bypanel') axis(2,cex.axis=binyaxcex)
					if(prs==nprs) axis(1,cex.axis=binxaxcex)
				} else {
					plot(1:3,1:3,type='n',axes=F,xlab='',ylab='')
				}
			} # if(any
			if(lat==nlat){ mtext(paste(prsmax[prs],' - ',prsmin[prs],sep=''),4,0.5,cex=0.75*binprscex) }
			if(prs==1){ mtext(paste(latmin[lat],' to ',latmax[lat],sep=''),3,0.5,cex=0.75*binlatcex) }
		} # loop on lat
	} # loop on prs
	mtext('Year Fraction',1,2,outer=T,cex=0.75)
	mtext(paste(species,' (',units,')',sep=''),2,2,outer=T,cex=0.75)
	mtext('Latitude (deg. N)',3,2,outer=T,cex=0.75)
	mtext('Pressure (hPa)',4,2,outer=T,cex=0.75)

	dev.off()


	# write out results
	print(Sys.time());print('outputting')
	write(latvals,paste(outdir,'/Seasonal_Amplitude_',species,flag,'.txt',sep=''),ncol=nlat)
	write(prsvals,paste(outdir,'/Seasonal_Amplitude_',species,flag,'.txt',sep=''),ncol=nprs,append=T)
	write(t(seasamp),paste(outdir,'/Seasonal_Amplitude_',species,flag,'.txt',sep=''),ncol=nprs,append=T)

	write(latvals,paste(outdir,'/Annual_Mean_',species,flag,'.txt',sep=''),ncol=nlat)
	write(prsvals,paste(outdir,'/Annual_Mean_',species,flag,'.txt',sep=''),ncol=nprs,append=T)
	write(t(annmean),paste(outdir,'/Annual_Mean_',species,flag,'.txt',sep=''),ncol=nprs,append=T)

	write(latvals,paste(outdir,'/Seasonal_Phase_',species,flag,'.txt',sep=''),ncol=nlat)
	write(prsvals,paste(outdir,'/Seasonal_Phase_',species,flag,'.txt',sep=''),ncol=nprs,append=T)
	write(t(seasphs),paste(outdir,'/Seasonal_Phase_',species,flag,'.txt',sep=''),ncol=nprs,append=T)

	write(latvals,paste(outdir,'/Count_',species,flag,'.txt',sep=''),ncol=nlat)
	write(prsvals,paste(outdir,'/Count_',species,flag,'.txt',sep=''),ncol=nprs,append=T)
	write(t(count),paste(outdir,'/Count_',species,flag,'.txt',sep=''),ncol=nprs,append=T)

	# calculating column means
	print(Sys.time());print('calculating column means')
	# this method averages previously binned (e.g. 100 hPa) averages for each column
	if(aggby=='day'){
		colval=aggregate(binval$x[prsvals[binval$prs]>coltop],by=list(lat=binval$lat[prsvals[binval$prs]>coltop],yday=binval$yday[prsvals[binval$prs]>coltop]),mean,na.rm=T)
		colwt=aggregate(!is.na(binval$x[prsvals[binval$prs]>coltop]),by=list(lat=binval$lat[prsvals[binval$prs]>coltop],yday=binval$yday[prsvals[binval$prs]>coltop]),sum)
	} else if(aggby=='campflt'){
		colval=aggregate(binval$x[prsvals[binval$prs]>coltop],by=list(lat=binval$lat[prsvals[binval$prs]>coltop],camp=binval$camp[prsvals[binval$prs]>coltop],flt=binval$flt[prsvals[binval$prs]>coltop]),mean,na.rm=T)
		colwt=aggregate(!is.na(binval$x[prsvals[binval$prs]>coltop]),by=list(lat=binval$lat[prsvals[binval$prs]>coltop],camp=binval$camp[prsvals[binval$prs]>coltop],flt=binval$flt[prsvals[binval$prs]>coltop]),sum)
		colval$yday=aggregate(binval$yday[prsvals[binval$prs]>coltop],by=list(lat=binval$lat[prsvals[binval$prs]>coltop],camp=binval$camp[prsvals[binval$prs]>coltop],flt=binval$flt[prsvals[binval$prs]>coltop]),mean,na.rm=T)$x
	} else if(aggby=='campslc'){
		colval=aggregate(binval$x[prsvals[binval$prs]>coltop],by=list(lat=binval$lat[prsvals[binval$prs]>coltop],camp=binval$camp[prsvals[binval$prs]>coltop],slc=binval$slc[prsvals[binval$prs]>coltop]),mean,na.rm=T)
		colwt=aggregate(!is.na(binval$x[prsvals[binval$prs]>coltop]),by=list(lat=binval$lat[prsvals[binval$prs]>coltop],camp=binval$camp[prsvals[binval$prs]>coltop],slc=binval$slc[prsvals[binval$prs]>coltop]),sum)
		colval$yday=aggregate(binval$yday[prsvals[binval$prs]>coltop],by=list(lat=binval$lat[prsvals[binval$prs]>coltop],camp=binval$camp[prsvals[binval$prs]>coltop],slc=binval$slc[prsvals[binval$prs]>coltop]),mean,na.rm=T)$x
	}
	print(paste('min,med,max colwt =',paste(round(c(min(colwt$x/(1000-coltop)*prsspan,na.rm=T),median(colwt$x/(1000-coltop)*prsspan,na.rm=T),max(colwt$x/(1000-coltop)*prsspan,na.rm=T)),2),collapse=', ')))
	## colwt is number of pressure bins filled on a particular day/flight/slice and colwtco is intended to remove partial columns

	print(paste('removing',sum(colwt$x/((1000-coltop)/prsspan)<colwtco),'of',nrow(colval),'column means'))
	colval$x[colwt$x/((1000-coltop)/prsspan)<colwtco]=NA

	write(colnames(colval),paste(outdir,'/Col_Values_',species,flag,'.txt',sep=''),ncol=ncol(colval))
	write(rbind(latvals[colval$lat],t(colval[,2:ncol(colval)])),paste(outdir,'/Col_Values_',species,flag,'.txt',sep=''),append=T,ncol=ncol(colval))

	## alternate method averages daily bin fits
	# used here to set ylim, but also checked and was identical to averaging bin averages
	if(sum(prsvals>coltop)>1){
		coldaily=apply(alldaily[,prsvals>coltop,],c(1,3),'mean',na.rm=T)
	} else {
		coldaily=alldaily
	}


	pngw=1800;pngh=1200;psize=24
	png(paste(outdir,'/Col_Fits_',species,flag,'.png',sep=''),width=pngw,height=pngh/5,pointsize=psize)
	par(mfrow=c(1,nlat))
	par(mar=c(0.1,0.1,0.1,0.1))
	par(oma=c(4,4,4,4))
	par(mgp=c(2.5,0.75,0))

	print(Sys.time());print('fitting column harmonics')
	# generage harmonics
	yrfrac=(colval$yday+0.5)/365*2*pi # could use difftime to do more precisely
	harm=cbind(cos(yrfrac),sin(yrfrac),cos(2*yrfrac),sin(2*yrfrac),cos(3*yrfrac),sin(3*yrfrac),cos(4*yrfrac),sin(4*yrfrac))[,1:(nharm*2)]

	seasamp=rep(NA,nlat)
	annmean=rep(NA,nlat)
	seasphs=rep(NA,nlat)

	if(ylim1=='fixed'|ylim1=='bylat') ylm=range(c(coldaily,colval$x),na.rm=T)

	for(lat in c(1:nlat)){
		numpts=sum(!is.na(colval$x[colval$lat==lat]))
		numcamp=length(unique(colval$camp[colval$lat==lat&!is.na(colval$x)]))

		if(numpts>nharm*2&(!reqncamp|numcamp>nharm*2)){
			if(ylim1=='bypanel') ylm=range(colval$x[colval$lat==lat],na.rm=T)
			# fit harmonics
			coefs=lm(colval$x[colval$lat==lat] ~ harm[colval$lat==lat,])$coef
			# regenerate harmonics
			daily=drop(coefs%*%t(hdays))
			seasamp[lat]=diff(range(daily))
			annmean[lat]=coefs[1]
			seasphs[lat]=seq(1,365)[daily==min(daily)][1] # in case of ties, pick first one

			plot((colval$yday[colval$lat==lat]+0.5)/365,colval$x[colval$lat==lat],main='',xlab='',ylab='',axes=F,cex.main=0.75,cex.axis=0.25,cex.lab=1.00,ylim=ylm,xlim=c(0,1),type='n')
			lines(seq(0.5,364.5)/365,daily,col='black')
			# lines(seq(0.5,364.5)/365,coldaily[lat,],col='cyan') # checked and was identical
			points((colval$yday[colval$lat==lat&colval$camp>20]+0.5)/365,colval$x[colval$lat==lat&colval$camp>20],col='dark blue') # ATom
		} else {
			if(grepl('MED_m_',species)|grepl('PFP_m_',species)){ # for flask-insitu differences, just fit a constant
				colmean=mean(colval$x[colval$lat==lat],na.rm=T)
				daily=rep(colmean,365)
				plot((colval$yday[colval$lat==lat]+0.5)/365,colval$x[colval$lat==lat],main='',xlab='',ylab='',axes=F,cex.main=0.75,cex.axis=0.25,cex.lab=1.00,ylim=ylm,xlim=c(0,1),type='n')
				lines(seq(0.5,364.5)/365,daily,col='black')
				points((colval$yday[colval$lat==lat&colval$camp>20]+0.5)/365,colval$x[colval$lat==lat&colval$camp>20],col='dark blue') # ATom
				# all colval$camp values = 1  2  3  4  5 11 21 22 23 24
			} else {
				plot(1:3,1:3,type='n',axes=F,xlab='',ylab='')
			}
		} # if(any
		box()
		if(lat==1&ylim1!='bypanel') axis(2,cex.axis=1.0)
		axis(1,cex.axis=1.0)
		mtext(paste(latmin[lat],' to ',latmax[lat],sep=''),3,0.5,cex=0.75)
	} # loop on lat

	mtext('Year Fraction',1,2,outer=T,cex=0.75)
	mtext(paste(species,' (',units,')',sep=''),2,2,outer=T,cex=0.75)
	mtext('Latitude (deg. N)',3,2,outer=T,cex=0.75)

	dev.off()


	# calculate curtain average metrics

	latwt=cos(latvals*pi/180)
	massatm=5.1352E18 # (dry kg) Trenberth and Smith, J. Climate, 2005
	molesatm=massatm*1E3/28.97
	rearth=6.371E6 # (m)
	# area of sphere = 4*pi*rearth^2, area of zone = 2*pi*rearth*h where h from equator = rearth*sin(lat)
	nmoles=molesatm*2*pi*rearth*(rearth-rearth*sin(calatco*pi/180))/(4*pi*rearth^2)*(1013.25-coltop)/1013.25 ## need to do a better job of mapping to actual surface pressure 

	pngw=1200;pngh=900;psize=12
	png(paste(outdir,'/CA_Fit_',species,flag,'.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	par(mfrow=c(1,2))
	par(mgp=c(2.5,0.75,0))

	if(-1*calatco>latrange[1]&-1*calatco<latrange[2]){
		# trim south of -calatco (e.g. > 20 S)
		latwt[latvals>(-1*calatco)]=0 ## assumes even division between bins at calatco
		cadaily=drop((latwt%*%coldaily)/sum(latwt))
		caamp=max(cadaily)-min(cadaily)
		print(paste('Curtain average south of ',calatco,' S seasonal amplitude = ',round(caamp,3),' ',units,sep=''))
		plot(seq(0.5,364.5)/365,cadaily,main='Curtain Average South',xlab=paste(species,' (',units,')',sep=''),ylab='Year Fraction',cex.main=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(0,1))
	}


	if(calatco>latrange[1]&calatco<latrange[2]){
		latwt=cos(latvals*pi/180)
		# trim north of calatco (e.g. > 20 N)
		latwt[latvals<calatco]=0 ## assumes even division between bins at calatco
		cadaily=drop((latwt%*%coldaily)/sum(latwt))
		caamp=max(cadaily)-min(cadaily)
		print(paste('Curtain average north of ',calatco,' N seasonal amplitude = ',round(caamp,3),' ',units,sep=''))
		plot(seq(0.5,364.5)/365,cadaily,main='Curtain Average North',xlab=paste(species,' (',units,')',sep=''),ylab='Year Fraction',cex.main=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(0,1))
	}

	dev.off()

	if(substr(species,1,3)=='CO2'&calatco>latrange[1]&calatco<latrange[2]){
		cadailyPgC=cadaily*1E-6*nmoles*12E-15 # assumes ppm
		caampPgC=max(cadailyPgC)-min(cadailyPgC)
		print(paste('Curtain average north of ',calatco,' seasonal amplitude = ',round(caampPgC,3),' PgC',sep=''))
	}


	# write out results
	print(Sys.time());print('outputting')
	write(rbind(latvals,seasamp),paste(outdir,'/Seasonal_Amplitude_',species,flag,'_Column.txt',sep=''),ncol=2)
	write(rbind(latvals,annmean),paste(outdir,'/Annual_Mean_',species,flag,'_Column.txt',sep=''),ncol=2)
	write(rbind(latvals,seasphs),paste(outdir,'/Seasonal_Phase_',species,flag,'_Column.txt',sep=''),ncol=2)

	# calculating surface values
	print(Sys.time());print('calculating surface values')
	if(aggby=='day'){
		surfval=aggregate(binval$x[prsvals[binval$prs]>=surfco],by=list(lat=binval$lat[prsvals[binval$prs]>=surfco],yday=binval$yday[prsvals[binval$prs]>=surfco]),mean,na.rm=T)
	} else if(aggby=='campflt'){
		surfval=aggregate(binval$x[prsvals[binval$prs]>=surfco],by=list(lat=binval$lat[prsvals[binval$prs]>=surfco],camp=binval$camp[prsvals[binval$prs]>=surfco],flt=binval$flt[prsvals[binval$prs]>=surfco]),mean,na.rm=T)
		surfval$yday=aggregate(binval$yday[prsvals[binval$prs]>=surfco],by=list(lat=binval$lat[prsvals[binval$prs]>=surfco],camp=binval$camp[prsvals[binval$prs]>=surfco],flt=binval$flt[prsvals[binval$prs]>=surfco]),mean,na.rm=T)$x
	} else if(aggby=='campslc'){
		surfval=aggregate(binval$x[prsvals[binval$prs]>=surfco],by=list(lat=binval$lat[prsvals[binval$prs]>=surfco],camp=binval$camp[prsvals[binval$prs]>=surfco],slc=binval$flt[prsvals[binval$prs]>=surfco]),mean,na.rm=T)
		surfval$yday=aggregate(binval$yday[prsvals[binval$prs]>=surfco],by=list(lat=binval$lat[prsvals[binval$prs]>=surfco],camp=binval$camp[prsvals[binval$prs]>=surfco],slc=binval$flt[prsvals[binval$prs]>=surfco]),mean,na.rm=T)$x
	}
	print(paste(nrow(surfval),'surface values'))

	print(Sys.time());print('fitting surface harmonics')
	# generage harmonics
	yrfrac=(surfval$yday+0.5)/365*2*pi # could use difftime to do more precisely
	harm=cbind(cos(yrfrac),sin(yrfrac),cos(2*yrfrac),sin(2*yrfrac),cos(3*yrfrac),sin(3*yrfrac),cos(4*yrfrac),sin(4*yrfrac))[,1:(nharm*2)]

	seasamp=rep(NA,nlat)
	annmean=rep(NA,nlat)
	seasphs=rep(NA,nlat)

	for(lat in c(1:nlat)){
		numpts=sum(!is.na(surfval$x[surfval$lat==lat]))
		numcamp=length(unique(surfval$camp[surfval$lat==lat&!is.na(surfval$x)]))
		if(numpts>nharm*2&(!reqncamp|numcamp>nharm*2)){
			# fit harmonics
			coefs=lm(surfval$x[surfval$lat==lat] ~ harm[surfval$lat==lat,])$coef
			# regenerate harmonics
			daily=drop(coefs%*%t(hdays))
			seasamp[lat]=diff(range(daily))
			annmean[lat]=coefs[1]
			seasphs[lat]=seq(1,365)[daily==min(daily)][1] # in case of ties, pick first one
		} else if(grepl('MED_m_',species)|grepl('PFP_m_',species)){
			surfmean=mean(surfval$x[surfval$lat==lat],na.rm=T)
			daily=rep(surfmean,365)
			seasamp[lat]=0
			annmean[lat]=surfmean
			seasphs[lat]=1
		}
	} # loop on lat


	# write out results
	print(Sys.time());print('outputting')
	write(rbind(latvals,seasamp),paste(outdir,'/Seasonal_Amplitude_',species,flag,'_Surface.txt',sep=''),ncol=2)
	write(rbind(latvals,annmean),paste(outdir,'/Annual_Mean_',species,flag,'_Surface.txt',sep=''),ncol=2)
	write(rbind(latvals,seasphs),paste(outdir,'/Seasonal_Phase_',species,flag,'_Surface.txt',sep=''),ncol=2)

} # end of harmf function
