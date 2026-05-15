## program to read in flux_proc.r output and make all-experiment plot of flux vs flux (e.g. O vs L, NET vs T)


fvf=function(start=ISOdate(2015,1,1,tz='UTC'),end=ISOdate(2020,12,31,tz='UTC')){

	library('ncdf4')
	library('RColorBrewer')
	cols=brewer.pal(5,'Set1')[c(1,3,2,4,5)]

	## Options:
	lonlim=c(-180,180) # for flux calculations
	exps=c('IS','OG','LNLG','LNLGOGIS','LNLGIS') # ,'unopt')
	mods=c('AMES','BAKER','CAMS','CMS-FLUX','COLA','CSU','CT','NIES','OU','TM5-4DVAR','UT','WEIR','WOMBAT') # all 14 except JHU
	modnames=c('Ames','Baker','CAMS','CMS-Flux','COLA','CSU','CT','NIES','OU','TM5-4DVAR','UT','LoFI','WOMBAT') # all 14 except JHU
	# JHU ocean fluxes are all identical and ~ 1 PgC/yr high

	#       Brad says: "IS is LoFI, which I consider a prior/“free-running”. It’s been categorized
	#       under IS because we use the NOAA MBL growth rate to adjust my global total fluxes. This is
	#       what’s described in our recent ACP paper: http://dx.doi.org/10.5194/acp-21-9609-2021
	#       LNLGOGIS is our OCO-2 assimilation, which we call OCO-2/GEOS. It uses LoFI as a prior.
	#       This is described in our recent Science Advances paper:
	#       http://dx.doi.org/10.1126/sciadv.abf9415 "
	#       But for ATom file, all experiments for LoFI are the same (LoFI) and all experiments for
	#       M2CC are the same (OCO-2/GEOS)

	starttxt=strftime(start,format='%Y-%m-%d')
	endtxt=strftime(end,format='%Y-%m-%d')


	cacalclist=c( # slat, nlat, limtxt, plot?
	-90,-80,'90S-80S',
	-80,-70,'80S-70S',
	-70,-60,'70S-60S',
	-60,-50,'60S-50S',
	-50,-40,'50S-40S',
	-40,-30,'40S-30S',
	-30,-20,'30S-20S',
	-20,-10,'20S-10S',
	-10,0,'10S-EQ',
	0,10,'EQ-10N',
	10,20,'10N-20N',
	20,30,'20N-30N',
	30,40,'30N-40N',
	40,50,'40N-50N',
	50,60,'50N-60N',
	60,70,'60N-70N',
	70,80,'70N-80N',
	80,90,'80N-90N'
	)
	# global means are calculated by summing these

	caspecs=matrix(cacalclist,byrow=T,ncol=3)

	cols=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2')) # 17 diff colors
	cols[6]='gold' # yellow too hard to see
	library('colorspace')
	bgcols=lighten(cols,0.75)

	lats=seq(-85,85,10)

	# calculate AGR
	noaatrend=read.csv('../LOCALDATA/co2_mm_gl.csv',skip=40) # from: https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_gl.csv
	colnames(noaatrend)=c('year','month','decimal','average','average_unc','trend','trend_unc')
	noaadt=ISOdate(noaatrend$year,noaatrend$month,15,tz='UTC')
	startco2=approx(noaadt,noaatrend$trend,start)$y
	endco2=approx(noaadt,noaatrend$trend,end+86400)$y
	years=as.numeric(difftime(end+86400,start,units='days'))/365.25 # not allowing exactly for leap years
	agr=(endco2-startco2)*2.124/years
	startunc=approx(noaadt,noaatrend$trend_unc,start)$y ## equivalent to GCB average of Dec and Jan
	endunc=approx(noaadt,noaatrend$trend_unc,end+86400)$y
	agrunc=((startunc^2+endunc^2)^0.5)*2.124/years
	# agrunc = 0.03 (for 6 years)
	# checked for consistency with GCB (0.02 for 10 years)


	## loop on caspecs, calculating global and zonal means for 6/2016-5/2018

	write('exp mod ocean land fossil',paste('oco2_v10mip_global_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''))
	write('exp mod ocean land fossil',paste('oco2_v10mip_net_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''))
	write('exp mod ocean land fossil',paste('oco2_v10mip_trop_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''))
	write('exp mod ocean land fossil',paste('oco2_v10mip_set_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''))
	write('exp mod ocean land fossil',paste('oco2_v10mip_et_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''))
	latbounds=c(-20,20)

	# loop on experiment
	globalmeanfluxes=NULL
	netmeanfluxes=NULL
	tropmeanfluxes=NULL
	setmeanfluxes=NULL
	etmeanfluxes=NULL

	for(exp in exps){
		print(exp)

		# loop on models
		for(mod in mods){
			print(mod)

			if(mod!='WEIR'|(exp=='IS'&exp!='unopt')){ # WEIR fluxes are the same for all experiments, but only apply to LoFI simulation which is considered IS (but no prior)

				globalmeanocean=0; globalmeanland=0; globalmeanfoss=0
				netmeanocean=0; netmeanland=0; netmeanfoss=0
				tropmeanocean=0; tropmeanland=0; tropmeanfoss=0
				setmeanocean=0; setmeanland=0; setmeanfoss=0
				etmeanocean=0; etmeanland=0; etmeanfoss=0
				for(i in c(1:nrow(caspecs))){
					limtxt=caspecs[i,3]
					if(file.exists(paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''))){

							monflux=read.table(paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''),header=T)
							monfluxdt=as.POSIXlt(ISOdate(monflux$year,monflux$month,15,0),tz='UTC') # midnight on 15th of each month
							globalmeanocean=globalmeanocean+mean(monflux$ocean[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
							globalmeanland=globalmeanland+mean(monflux$land[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
							globalmeanfoss=globalmeanfoss+mean(monflux$foss[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
							if(lats[i]>latbounds[2]){ # NET
								netmeanocean=netmeanocean+mean(monflux$ocean[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
								netmeanland=netmeanland+mean(monflux$land[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
								netmeanfoss=netmeanfoss+mean(monflux$foss[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
							} 
							if(lats[i]<latbounds[1]){ # SET
								setmeanocean=setmeanocean+mean(monflux$ocean[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
								setmeanland=setmeanland+mean(monflux$land[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
								setmeanfoss=setmeanfoss+mean(monflux$foss[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
							} 
							if(lats[i]>latbounds[2]|lats[i]<latbounds[1]){ # ET
								etmeanocean=etmeanocean+mean(monflux$ocean[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
								etmeanland=etmeanland+mean(monflux$land[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
								etmeanfoss=etmeanfoss+mean(monflux$foss[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
							}
							if(lats[i]<latbounds[2]&lats[i]>latbounds[1]){ # Trop
								tropmeanocean=tropmeanocean+mean(monflux$ocean[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
								tropmeanland=tropmeanland+mean(monflux$land[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
								tropmeanfoss=tropmeanfoss+mean(monflux$foss[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
							}

					}

				} # loop on caspecs

				globalmeanfluxes=rbind(globalmeanfluxes,c(exp,mod,round(globalmeanocean,3),round(globalmeanland,3),round(globalmeanfoss,3)))
				netmeanfluxes=rbind(netmeanfluxes,c(exp,mod,round(netmeanocean,3),round(netmeanland,3),round(netmeanfoss,3)))
				tropmeanfluxes=rbind(tropmeanfluxes,c(exp,mod,round(tropmeanocean,3),round(tropmeanland,3),round(tropmeanfoss,3)))
				setmeanfluxes=rbind(setmeanfluxes,c(exp,mod,round(setmeanocean,3),round(setmeanland,3),round(setmeanfoss,3)))
				etmeanfluxes=rbind(etmeanfluxes,c(exp,mod,round(etmeanocean,3),round(etmeanland,3),round(etmeanfoss,3)))

			} # only use LoFI IS

		} # loop on model

	} # loop on experiment

	write(t(globalmeanfluxes),paste('oco2_v10mip_global_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),ncol=5,append=T)
	write(t(netmeanfluxes),paste('oco2_v10mip_net_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),ncol=5,append=T)
	write(t(tropmeanfluxes),paste('oco2_v10mip_trop_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),ncol=5,append=T)
	write(t(setmeanfluxes),paste('oco2_v10mip_set_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),ncol=5,append=T)
	write(t(etmeanfluxes),paste('oco2_v10mip_et_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),ncol=5,append=T)

	## do also for global unopt
	exp='unopt'
	priormeanfluxes=NULL

	# loop on models
	for(mod in mods){

		if(mod!='WEIR'){ # no prior

			priormeanocean=0
			priormeanland=0
			priormeanfoss=0
			for(i in c(1:nrow(caspecs))){
				limtxt=caspecs[i,3]
				if(file.exists(paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''))){
						monflux=read.table(paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''),header=T)
						monfluxdt=as.POSIXlt(ISOdate(monflux$year,monflux$month,15,0),tz='UTC') # midnight on 15th of each month
						priormeanocean=priormeanocean+mean(monflux$ocean[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
						priormeanland=priormeanland+mean(monflux$land[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
						priormeanfoss=priormeanfoss+mean(monflux$foss[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0])
				}
			} # loop on caspecs
			priormeanfluxes=rbind(priormeanfluxes,c(exp,mod,round(priormeanocean,3),round(priormeanland,3),round(priormeanfoss,3)))

		} # no WEIR prior

	} # loop on model
	write('exp mod ocean land fossil',paste('oco2_v10mip_prior_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''))
	write(t(priormeanfluxes),paste('oco2_v10mip_prior_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),ncol=5,append=T)


	# read back in:
	globalmeanfluxes=read.table(paste('oco2_v10mip_global_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),header=T,stringsAsFactors=F)
	netmeanfluxes=read.table(paste('oco2_v10mip_net_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),header=T,stringsAsFactors=F)
	tropmeanfluxes=read.table(paste('oco2_v10mip_trop_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),header=T,stringsAsFactors=F)
	setmeanfluxes=read.table(paste('oco2_v10mip_set_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),header=T,stringsAsFactors=F)
	etmeanfluxes=read.table(paste('oco2_v10mip_et_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),header=T,stringsAsFactors=F)
	priormeanfluxes=read.table(paste('oco2_v10mip_prior_fluxes_',starttxt,'_to_',endtxt,'.txt',sep=''),header=T,stringsAsFactors=F)

	globalmeanfluxes$total=globalmeanfluxes$ocean+globalmeanfluxes$land+globalmeanfluxes$foss
	netmeanfluxes$total=netmeanfluxes$ocean+netmeanfluxes$land+netmeanfluxes$foss
	tropmeanfluxes$total=tropmeanfluxes$ocean+tropmeanfluxes$land+tropmeanfluxes$foss
	setmeanfluxes$total=setmeanfluxes$ocean+setmeanfluxes$land+setmeanfluxes$foss
	etmeanfluxes$total=etmeanfluxes$ocean+etmeanfluxes$land+etmeanfluxes$foss
	priormeanfluxes$total=priormeanfluxes$ocean+priormeanfluxes$land+priormeanfluxes$foss

	# calculate FF, uncertainty and AGR+FF combined uncertainty
	foss=globalmeanfluxes$foss[1] # same for all models and experiments ### not true for other MIPs
	fossunc=0.05*foss # 5% is from GCB 2024 (Andres, 2012)
	agrfossunc=(agrunc^2+fossunc^2)^0.5 
	# 0.487 without AGR, 0.492 with

	expcols=c(brewer.pal(5,'Set1'),'black')

	png(paste('oco_v10mip_fluxes_ocean_vs_land_byExp_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)
	par(mar=c(5, 5, 4, 1)+0.1)

	plot(globalmeanfluxes$land,globalmeanfluxes$ocean,xlab=expression(paste('Global Land Flux (PgC ',yr^-1,')')),ylab=expression(paste('Global Ocean Flux (PgC ',yr^-1,')')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('a. Global Ocean vs. Land Fluxes'),type='n')
	xcorners=c(-5,1,1,-5,-5)
	ycorners=-1*xcorners+agr-foss+c(-1*agrfossunc,-1*agrfossunc,agrfossunc,agrfossunc,-1*agrfossunc)
	polygon(xcorners,ycorners,col='gray75',border=NA)
	abline(agr-foss,-1,lwd=3)
	for(mod in mods){
		points(globalmeanfluxes$land[globalmeanfluxes$mod==mod],globalmeanfluxes$ocean[globalmeanfluxes$mod==mod],col=expcols[match(globalmeanfluxes$exp[globalmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white')
	}

	for(exp in exps){
		expmeanocean=mean(globalmeanfluxes$ocean[globalmeanfluxes$exp==exp])
		expsdocean=sd(globalmeanfluxes$ocean[globalmeanfluxes$exp==exp])
		expmeanland=mean(globalmeanfluxes$land[globalmeanfluxes$exp==exp])
		expsdland=sd(globalmeanfluxes$land[globalmeanfluxes$exp==exp])
		points(expmeanland,expmeanocean,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeanland-expsdland,expmeanocean,expmeanland+expsdland,expmeanocean,lwd=4,col=expcols[match(exp,exps)])
		segments(expmeanland,expmeanocean-expsdocean,expmeanland,expmeanocean+expsdocean,lwd=4,col=expcols[match(exp,exps)])
	}
	legend('topright',modnames,col='black',pt.lwd=3,pch=c(1:length(mods)),pt.bg='white',bty='n',title='Model:',title.col='black',title.adj=0.1)
	legend('bottomleft',exps,text.col=cols[c(1:length(exps))],pch=16,col=cols[c(1:length(exps))],pt.cex=2.0,pt.lwd=3,pt.bg='white',lty=NA,bty='n',title='Experiment:',title.col='black',title.adj=0.1,lwd=1) # crashes without lwd=1
	legend('top',legend=substitute(AGR - FF == a %+-% b,list(a=round(agr-foss,2),b=round(agrfossunc,2))),pch=NA,lty=1,bty='n',lwd=3)
	dev.off()

	## EPS Version:

	# --- Dimensions and Scaling Logic ---
	width_cm <- 8.9
	width_in <- width_cm / 2.54
	sc <- 0.21 # Scaling factor for absolute line widths
	new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)

	file_name <- paste('oco_v10mip_fluxes_ocean_vs_land_byExp_', starttxt, '_to_', endtxt, '.eps', sep='')

	cairo_ps(filename = file_name,
		 width = width_in,
		 height = width_in,
		 pointsize = new_pointsize)

	# --- Plotting Code ---
	par(mar=c(5, 5, 4, 1) + 0.1)

	# cex.axis/lab/main are multipliers of pointsize; 1.5 is the original ratio
	plot(globalmeanfluxes$land, globalmeanfluxes$ocean,
	     xlab=expression(paste('Global Land Flux (PgC ', yr^-1, ')')),
	     ylab=expression(paste('Global Ocean Flux (PgC ', yr^-1, ')')),
	     cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
	     main=expression('a. Global Ocean vs. Land Fluxes'), type='n',axes=F)

	xcorners = c(-5, 1, 1, -5, -5)
	ycorners = -1*xcorners + agr - foss + c(-1*agrfossunc, -1*agrfossunc, agrfossunc, agrfossunc, -1*agrfossunc)

	polygon(xcorners, ycorners, col='gray75', border=NA)
	abline(agr-foss, -1, lwd=3 * sc) # Keep lwd scaled

	box(lwd=0.5)
	axis(1,cex.axis=1.5,lwd=0.5)
	axis(2,cex.axis=1.5,lwd=0.5)

	for(mod in mods){

	    points(globalmeanfluxes$land[globalmeanfluxes$mod==mod],
		   globalmeanfluxes$ocean[globalmeanfluxes$mod==mod],
		   col=expcols[match(globalmeanfluxes$exp[globalmeanfluxes$mod==mod],exps)],
		   cex=2,
		   lwd=3 * sc,
		   pch=match(mod,mods),
		   bg='white') 

	}

	for(exp in exps){
		expmeanocean=mean(globalmeanfluxes$ocean[globalmeanfluxes$exp==exp])
		expsdocean=sd(globalmeanfluxes$ocean[globalmeanfluxes$exp==exp])
		expmeanland=mean(globalmeanfluxes$land[globalmeanfluxes$exp==exp])
		expsdland=sd(globalmeanfluxes$land[globalmeanfluxes$exp==exp])
		points(expmeanland,expmeanocean,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeanland-expsdland,expmeanocean,expmeanland+expsdland,expmeanocean,lwd=4*sc,col=expcols[match(exp,exps)])
		segments(expmeanland,expmeanocean-expsdocean,expmeanland,expmeanocean+expsdocean,lwd=4*sc,col=expcols[match(exp,exps)])
	}

	legend('topright',modnames,col='black',pt.lwd=3*sc,pch=c(1:length(mods)),pt.bg='white',bty='n',title='Model:',title.col='black',title.adj=0.1)

	legend('bottomleft',exps,text.col=cols[c(1:length(exps))],pch=16,col=cols[c(1:length(exps))],pt.cex=2.0,pt.lwd=3*sc,pt.bg='white',lty=NA,bty='n',title='Experiment:',title.col='black',title.adj=0.1,lwd=1)

	legend('top',legend=substitute(AGR - FF == a %+-% b,list(a=round(agr-foss,2),b=round(agrfossunc,2))),pch=NA,lty=1,bty='n',lwd=3*sc)

	dev.off()


	modmeanocean=NULL
	modsdocean=NULL
	modrngocean=NULL
	modmeanland=NULL
	modsdland=NULL
	modrngland=NULL
	for(mod in mods){
		modmeanocean=c(modmeanocean,mean(globalmeanfluxes$ocean[globalmeanfluxes$mod==mod]))
		modsdocean=c(modsdocean,sd(globalmeanfluxes$ocean[globalmeanfluxes$mod==mod]))
		modrngocean=c(modrngocean,diff(range(globalmeanfluxes$ocean[globalmeanfluxes$mod==mod],na.rm=T)))
		modmeanland=c(modmeanland,mean(globalmeanfluxes$land[globalmeanfluxes$mod==mod]))
		modsdland=c(modsdland,sd(globalmeanfluxes$land[globalmeanfluxes$mod==mod]))
		modrngland=c(modrngland,diff(range(globalmeanfluxes$land[globalmeanfluxes$mod==mod],na.rm=T)))
	}
	modsds=cbind(mods,round(modsdocean,2),round(modsdland,2))
	modsds=modsds[!is.na(modsdocean),]
	modsds=data.frame(modsds)
	colnames(modsds)=c('mod','ocean','land')
	print(modsds)

	modrngs=cbind(mods,round(modrngocean,2),round(modrngland,2))
	modrngs=modrngs[!is.na(modrngocean)&!is.infinite(modrngocean)&modrngocean!=0,]
	modrngs=data.frame(modrngs)
	colnames(modrngs)=c('mod','ocean','land')
	print(modrngs)

	png(paste('oco_v10mip_fluxes_ocean_vs_land_byMod_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)
	par(mar=c(5, 5, 4, 1)+0.1)

	plot(globalmeanfluxes$land,globalmeanfluxes$ocean,xlab=expression(paste('Global Land Flux (PgC ',yr^-1,')')),ylab=expression(paste('Global Ocean Flux (PgC ',yr^-1,')')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('a. Global Ocean vs. Land Fluxes'),type='n')
	modord=order(modsdland,decreasing=T)
	xcorners=c(-5,1,1,-5,-5)
	ycorners=-1*xcorners+agr-foss+c(-1*agrfossunc,-1*agrfossunc,agrfossunc,agrfossunc,-1*agrfossunc)
	polygon(xcorners,ycorners,col='gray75',border=NA)
	abline(agr-foss,-1,lwd=3)
	for(mod in mods[modord]){
		expord=order(globalmeanfluxes$land[globalmeanfluxes$mod==mod])
		exppch=match(globalmeanfluxes$exp[globalmeanfluxes$mod==mod][expord],exps)+20
		exppch[exppch==26]=21 # added for testing unopt
		points(globalmeanfluxes$land[globalmeanfluxes$mod==mod][expord],globalmeanfluxes$ocean[globalmeanfluxes$mod==mod][expord],col=cols[match(mod,mods)],cex=2,lwd=3,pch=exppch,bg=bgcols[match(mod,mods)],type='b') #,'white')
	}
	legend('topright',modnames,pch=NA,lty=NA,text.col=cols[c(1:length(mods))],bty='n',title='Model:',title.col='black',title.adj=0.5,lwd=1)
	legend('bottomleft',exps,col='black',pch=21:25,lty=NA,lwd=3,pt.cex=1.5,bty='n',pt.bg='white',title='Experiment:',title.col='black',title.adj=0.1)
	legend('top',legend=substitute(AGR - FF == a %+-% b,list(a=round(agr-foss,2),b=round(agrfossunc,2))),pch=NA,lty=1,bty='n',lwd=3)

	dev.off()

	## EPS Version:

	# --- Dimensions and Scaling Logic ---
	width_cm <- 8.9
	width_in <- width_cm / 2.54
	sc <- 0.21 # Scaling factor for absolute line widths
	new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)

	file_name <- paste('oco_v10mip_fluxes_ocean_vs_land_byMod_', starttxt, '_to_', endtxt, '.eps', sep='')

	cairo_ps(filename = file_name,
		 width = width_in,
		 height = width_in,
		 pointsize = new_pointsize)

	# --- Plotting Code ---
	par(mar=c(5, 5, 4, 1) + 0.1)

	# cex.axis/lab/main are multipliers of pointsize; 1.5 is the original ratio
	plot(globalmeanfluxes$land, globalmeanfluxes$ocean,
	     xlab=expression(paste('Global Land Flux (PgC ', yr^-1, ')')),
	     ylab=expression(paste('Global Ocean Flux (PgC ', yr^-1, ')')),
	     cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
	     main=expression('a. Global Ocean vs. Land Fluxes'), type='n',axes=F)

	modord = order(modsdland, decreasing=T)
	xcorners = c(-5, 1, 1, -5, -5)
	ycorners = -1*xcorners + agr - foss + c(-1*agrfossunc, -1*agrfossunc, agrfossunc, agrfossunc, -1*agrfossunc)

	polygon(xcorners, ycorners, col='gray75', border=NA)
	abline(agr-foss, -1, lwd=3 * sc) # Keep lwd scaled

	box(lwd=0.5)
	axis(1,cex.axis=1.5,lwd=0.5)
	axis(2,cex.axis=1.5,lwd=0.5)

	for(mod in mods[modord]){
	    expord = order(globalmeanfluxes$land[globalmeanfluxes$mod==mod])
	    exppch = match(globalmeanfluxes$exp[globalmeanfluxes$mod==mod][expord], exps) + 20
	    exppch[exppch==26] = 21

	    points(globalmeanfluxes$land[globalmeanfluxes$mod==mod][expord],
		   globalmeanfluxes$ocean[globalmeanfluxes$mod==mod][expord],
		   col=cols[match(mod, mods)],
		   cex=2,
		   lwd=3 * sc,
		   pch=exppch,
		   bg=bgcols[match(mod, mods)], type='b')
	}

	legend('topright', modnames, pch=NA, lty=NA, text.col=cols[c(1:length(mods))],
	       bty='n', title='Model:', title.col='black', title.adj=0.5, lwd=1)

	legend('bottomleft', exps, col='black', pch=21:25, lty=NA,
	       lwd=3 * sc, pt.cex=1.5,
	       bty='n', pt.bg='white', title='Experiment:', title.col='black', title.adj=0.1)

	legend('top', legend=substitute(AGR - FF == a %+-% b, list(a=round(agr-foss, 2), b=round(agrfossunc, 2))),
	       pch=NA, lty=1, bty='n', lwd=3 * sc)

	dev.off()


	png(paste('oco_v10mip_fluxes_total_vs_Exp_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24) # like Gaubert et al., 2019 Fig. 4b

	plot(match(globalmeanfluxes$exp,exps),globalmeanfluxes$total,xlab='Experiment',ylab='Global Total Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Global Total vs. Experiment',type='n')
	for(mod in mods){
		points(match(globalmeanfluxes$exp,exps)[globalmeanfluxes$mod==mod],globalmeanfluxes$total[globalmeanfluxes$mod==mod],col=expcols[match(globalmeanfluxes$exp[globalmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white')
	}

	for(exp in exps){
		expmeantotal=mean(globalmeanfluxes$total[globalmeanfluxes$exp==exp])
		expsdtotal=sd(globalmeanfluxes$total[globalmeanfluxes$exp==exp])
		points(match(exp,exps),expmeantotal,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(match(exp,exps),expmeantotal-expsdtotal,match(exp,exps),expmeantotal+expsdtotal,lwd=4,col=expcols[match(exp,exps)])
	}
	mtext(paste('SD =',round(sd(globalmeanfluxes$total),2)),1,-1.5)
	dev.off()


	png(paste('oco_v10mip_fluxes_trop_vs_net_total_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)

	plot(netmeanfluxes$total,tropmeanfluxes$total,xlab='NET Total Flux (PgC/yr)',ylab='Trop Total Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Trop vs. NET Total Fluxes',type='n')
	for(mod in mods){
		points(netmeanfluxes$total[netmeanfluxes$mod==mod],tropmeanfluxes$total[tropmeanfluxes$mod==mod],col=expcols[match(netmeanfluxes$exp[netmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white')
	}

	for(exp in exps){
		expmeantrop=mean(tropmeanfluxes$total[tropmeanfluxes$exp==exp])
		expsdtrop=sd(tropmeanfluxes$total[tropmeanfluxes$exp==exp])
		expmeannet=mean(netmeanfluxes$total[netmeanfluxes$exp==exp])
		expsdnet=sd(netmeanfluxes$total[netmeanfluxes$exp==exp])
		expmeanset=mean(setmeanfluxes$total[setmeanfluxes$exp==exp]) # not plotted
		expsdset=sd(setmeanfluxes$total[setmeanfluxes$exp==exp]) # not plotted
		points(expmeannet,expmeantrop,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeannet-expsdnet,expmeantrop,expmeannet+expsdnet,expmeantrop,lwd=4,col=expcols[match(exp,exps)])
		segments(expmeannet,expmeantrop-expsdtrop,expmeannet,expmeantrop+expsdtrop,lwd=4,col=expcols[match(exp,exps)])
		print(paste(exp,round(expmeannet,2),round(expsdnet/10^0.5,2),round(expmeantrop,2),round(expsdtrop/10^0.5,2),round(expmeanset,2),round(expsdset/10^0.5,2)))
	}

	dev.off()


	png(paste('oco_v10mip_fluxes_set_vs_net_total_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)

	plot(netmeanfluxes$total,setmeanfluxes$total,xlab='NET Total Flux (PgC/yr)',ylab='Trop Total Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP SET vs. NET Total Fluxes',type='n')
	for(mod in mods){
		points(netmeanfluxes$total[netmeanfluxes$mod==mod],setmeanfluxes$total[setmeanfluxes$mod==mod],col=expcols[match(netmeanfluxes$exp[netmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white')
	}
	print(paste('SET vs NET R^2 =',round(cor(netmeanfluxes$total,setmeanfluxes$total)^2,3)))
	print(paste('SET vs NET No OG R^2 =',round(cor(netmeanfluxes$total[netmeanfluxes$exp!='OG'&netmeanfluxes$exp!='LNLGOGIS'],setmeanfluxes$total[netmeanfluxes$exp!='OG'&netmeanfluxes$exp!='LNLGOGIS'])^2,3)))

	allexpmeansNET=NULL
	allexpmeansSET=NULL
	allexp=NULL
	for(exp in exps){
		expmeanset=mean(setmeanfluxes$total[setmeanfluxes$exp==exp])
		expsdset=sd(setmeanfluxes$total[setmeanfluxes$exp==exp])
		expmeannet=mean(netmeanfluxes$total[netmeanfluxes$exp==exp])
		expsdnet=sd(netmeanfluxes$total[netmeanfluxes$exp==exp])
		expmeanset=mean(setmeanfluxes$total[setmeanfluxes$exp==exp]) # not plotted
		expsdset=sd(setmeanfluxes$total[setmeanfluxes$exp==exp]) # not plotted
		points(expmeannet,expmeanset,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeannet-expsdnet,expmeanset,expmeannet+expsdnet,expmeanset,lwd=4,col=expcols[match(exp,exps)])
		segments(expmeannet,expmeanset-expsdset,expmeannet,expmeanset+expsdset,lwd=4,col=expcols[match(exp,exps)])
		print(paste(exp,round(expmeannet,2),round(expsdnet/10^0.5,2),round(expmeanset,2),round(expsdset/10^0.5,2),round(expmeanset,2),round(expsdset/10^0.5,2)))
		allexp=c(allexp,exp)
		allexpmeansNET=c(allexpmeansNET,expmeannet)
		allexpmeansSET=c(allexpmeansSET,expmeanset)
	}

	print(paste('SET vs NET Exp Means R^2 =',round(cor(allexpmeansNET,allexpmeansSET)^2,3)))
	print(paste('SET vs NET Exp Means No OG R^2 =',round(cor(allexpmeansNET[allexp!='OG'&allexp!='LNLGOGIS'],allexpmeansSET[allexp!='OG'&allexp!='LNLGOGIS'])^2,3)))

	dev.off()


	png(paste('oco_v10mip_fluxes_trop_vs_net_land_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)

	plot(netmeanfluxes$land,tropmeanfluxes$land,xlab='NET Land Flux (PgC/yr)',ylab='Trop Land Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Trop vs. NET Land Fluxes',type='n')
	for(mod in mods){
		points(netmeanfluxes$land[netmeanfluxes$mod==mod],tropmeanfluxes$land[tropmeanfluxes$mod==mod],col=expcols[match(netmeanfluxes$exp[netmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white')
	}

	for(exp in exps){
		expmeantrop=mean(tropmeanfluxes$land[tropmeanfluxes$exp==exp])
		expsdtrop=sd(tropmeanfluxes$land[tropmeanfluxes$exp==exp])
		expmeannet=mean(netmeanfluxes$land[netmeanfluxes$exp==exp])
		expsdnet=sd(netmeanfluxes$land[netmeanfluxes$exp==exp])
		points(expmeannet,expmeantrop,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeannet-expsdnet,expmeantrop,expmeannet+expsdnet,expmeantrop,lwd=4,col=expcols[match(exp,exps)])
		segments(expmeannet,expmeantrop-expsdtrop,expmeannet,expmeantrop+expsdtrop,lwd=4,col=expcols[match(exp,exps)])
	}

	dev.off()


	png(paste('oco_v10mip_fluxes_trop_vs_net_ocean_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)

	plot(netmeanfluxes$ocean,tropmeanfluxes$ocean,xlab='NET Ocean Flux (PgC/yr)',ylab='Trop Ocean Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Trop vs. NET Ocean Fluxes',type='n')
	for(mod in mods){
		points(netmeanfluxes$ocean[netmeanfluxes$mod==mod],tropmeanfluxes$ocean[tropmeanfluxes$mod==mod],col=expcols[match(netmeanfluxes$exp[netmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white')
	}

	for(exp in exps){
		expmeantrop=mean(tropmeanfluxes$ocean[tropmeanfluxes$exp==exp])
		expsdtrop=sd(tropmeanfluxes$ocean[tropmeanfluxes$exp==exp])
		expmeannet=mean(netmeanfluxes$ocean[netmeanfluxes$exp==exp])
		expsdnet=sd(netmeanfluxes$ocean[netmeanfluxes$exp==exp])
		points(expmeannet,expmeantrop,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeannet-expsdnet,expmeantrop,expmeannet+expsdnet,expmeantrop,lwd=4,col=expcols[match(exp,exps)])
		segments(expmeannet,expmeantrop-expsdtrop,expmeannet,expmeantrop+expsdtrop,lwd=4,col=expcols[match(exp,exps)])
	}

	dev.off()


	png(paste('oco_v10mip_fluxes_trop_vs_et_total_byExp_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)
	par(mar=c(5, 5, 4, 1)+0.1)

	plot(etmeanfluxes$total,tropmeanfluxes$total,xlab=expression(paste('Extratropical Total Flux (PgC ',yr^-1,')')),ylab=expression(paste('Tropical Total Flux (PgC ',yr^-1,')')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Tropical vs. Extratropical Total Fluxes'),type='n')
	xcorners=c(-1,7,7,-1,-1)
	ycorners=-1*xcorners+agr+c(-1*agrunc,-1*agrunc,agrunc,agrunc,-1*agrunc)
	polygon(xcorners,ycorners,col='gray75',border=NA)
	abline(agr,-1,lwd=3)
	for(mod in mods){
		points(etmeanfluxes$total[etmeanfluxes$mod==mod],tropmeanfluxes$total[tropmeanfluxes$mod==mod],col=expcols[match(etmeanfluxes$exp[etmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white') # expcols[match(etmeanfluxes$exp[etmeanfluxes$mod==mod],exps)]) # 'white')
	}

	for(exp in exps){
		expmeantrop=mean(tropmeanfluxes$total[tropmeanfluxes$exp==exp])
		expsdtrop=sd(tropmeanfluxes$total[tropmeanfluxes$exp==exp])
		expmeanet=mean(etmeanfluxes$total[etmeanfluxes$exp==exp])
		expsdet=sd(etmeanfluxes$total[etmeanfluxes$exp==exp])
		points(expmeanet,expmeantrop,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeanet-expsdet,expmeantrop,expmeanet+expsdet,expmeantrop,lwd=4,col=expcols[match(exp,exps)])
		segments(expmeanet,expmeantrop-expsdtrop,expmeanet,expmeantrop+expsdtrop,lwd=4,col=expcols[match(exp,exps)])
	}
	legend('topright',modnames,col='black',pt.lwd=3,pch=c(1:length(mods)),pt.bg='white',bty='n',title='Model:',title.col='black',title.adj=0.1)
	legend('bottomleft',exps,text.col=cols[c(1:length(exps))],pch=16,col=cols[c(1:length(exps))],pt.cex=2.0,pt.lwd=3,pt.bg='white',lty=NA,bty='n',title='Experiment:',title.col='black',title.adj=0.1,lwd=1)
	legend('top',legend=substitute(AGR == a %+-% b,list(a=round(agr,2),b=round(agrunc,2))),pch=NA,lty=1,bty='n',lwd=3)

	dev.off()

	## EPS version:

	# --- Dimensions and Scaling Logic ---
	width_cm <- 8.9
	width_in <- width_cm / 2.54
	sc <- 0.21 # Scaling factor for absolute line widths
	new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)

	file_name <- paste('oco_v10mip_fluxes_trop_vs_et_total_byExp_', starttxt, '_to_', endtxt, '.eps', sep='')

	cairo_ps(filename = file_name,
		 width = width_in,
		 height = width_in,
		 pointsize = new_pointsize)


	par(mar=c(5, 5, 4, 1)+0.1)

	plot(etmeanfluxes$total,tropmeanfluxes$total,xlab=expression(paste('Extratropical Total Flux (PgC ',yr^-1,')')),ylab=expression(paste('Tropical Total Flux (PgC ',yr^-1,')')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Tropical vs. Extratropical Total Fluxes'),type='n',axes=F)
	xcorners=c(-1,7,7,-1,-1)
	ycorners=-1*xcorners+agr+c(-1*agrunc,-1*agrunc,agrunc,agrunc,-1*agrunc)
	polygon(xcorners,ycorners,col='gray75',border=NA)
	abline(agr,-1,lwd=3*sc)

	box(lwd=0.5)
	axis(1,cex.axis=1.5,lwd=0.5)
	axis(2,cex.axis=1.5,lwd=0.5)

	for(mod in mods){
		points(etmeanfluxes$total[etmeanfluxes$mod==mod],tropmeanfluxes$total[tropmeanfluxes$mod==mod],col=expcols[match(etmeanfluxes$exp[etmeanfluxes$mod==mod],exps)],cex=2,lwd=3*sc,pch=match(mod,mods),bg='white') # expcols[match(etmeanfluxes$exp[etmeanfluxes$mod==mod],exps)]) # 'white')
	}

	for(exp in exps){
		expmeantrop=mean(tropmeanfluxes$total[tropmeanfluxes$exp==exp])
		expsdtrop=sd(tropmeanfluxes$total[tropmeanfluxes$exp==exp])
		expmeanet=mean(etmeanfluxes$total[etmeanfluxes$exp==exp])
		expsdet=sd(etmeanfluxes$total[etmeanfluxes$exp==exp])
		points(expmeanet,expmeantrop,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeanet-expsdet,expmeantrop,expmeanet+expsdet,expmeantrop,lwd=4*sc,col=expcols[match(exp,exps)])
		segments(expmeanet,expmeantrop-expsdtrop,expmeanet,expmeantrop+expsdtrop,lwd=4*sc,col=expcols[match(exp,exps)])
	}
	legend('topright',modnames,col='black',pt.lwd=3*sc,pch=c(1:length(mods)),pt.bg='white',bty='n',title='Model:',title.col='black',title.adj=0.1)
	legend('bottomleft',exps,text.col=cols[c(1:length(exps))],pch=16,col=cols[c(1:length(exps))],pt.cex=2.0,pt.lwd=3*sc,pt.bg='white',lty=NA,bty='n',title='Experiment:',title.col='black',title.adj=0.1,lwd=1)
	legend('top',legend=substitute(AGR == a %+-% b,list(a=round(agr,2),b=round(agrunc,2))),pch=NA,lty=1,bty='n',lwd=3*sc)

	dev.off()


	png(paste('oco_v10mip_fluxes_trop_vs_et_total_byMod_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)
	par(mar=c(5, 5, 4, 1)+0.1)

	modmeantrop=NULL
	modsdtrop=NULL
	modmeanet=NULL
	modsdet=NULL
	for(mod in mods){
		modmeantrop=c(modmeantrop,mean(tropmeanfluxes$total[tropmeanfluxes$mod==mod]))
		modsdtrop=c(modsdtrop,sd(tropmeanfluxes$total[tropmeanfluxes$mod==mod]))
		modmeanet=c(modmeanet,mean(etmeanfluxes$total[etmeanfluxes$mod==mod]))
		modsdet=c(modsdet,sd(etmeanfluxes$total[etmeanfluxes$mod==mod]))
	}

	plot(etmeanfluxes$total,tropmeanfluxes$total,xlab=expression(paste('Extratropical Total Flux (PgC ',yr^-1,')')),ylab=expression(paste('Tropical Total Flux (PgC ',yr^-1,')')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Tropical vs. Extratropical Total Fluxes'),type='n')
	xcorners=c(-1,7,7,-1,-1)
	ycorners=-1*xcorners+agr+c(-1*agrunc,-1*agrunc,agrunc,agrunc,-1*agrunc)
	polygon(xcorners,ycorners,col='gray75',border=NA)
	abline(agr,-1,lwd=3)
	modord=order(modsdet,decreasing=T)
	for(mod in mods[modord]){
		expord=order(etmeanfluxes$total[etmeanfluxes$mod==mod])
		exppch=match(etmeanfluxes$exp[etmeanfluxes$mod==mod][expord],exps)+20
		exppch[exppch==26]=21 # added for testing unopt
		points(etmeanfluxes$total[etmeanfluxes$mod==mod][expord],tropmeanfluxes$total[tropmeanfluxes$mod==mod][expord],col=cols[match(mod,mods)],cex=2,lwd=3,pch=exppch,bg=bgcols[match(mod,mods)],type='b') # ,'white')
	}
	legend('topright',modnames,pch=NA,lty=NA,text.col=cols[c(1:length(mods))],bty='n',title='Model:',title.col='black',title.adj=0.5,lwd=1)
	legend('bottomleft',exps,col='black',pch=21:25,lty=NA,lwd=3,pt.cex=1.5,bty='n',pt.bg='white',title='Experiment:',title.col='black',title.adj=0.1) 
	legend('top',legend=substitute(AGR == a %+-% b,list(a=round(agr,2),b=round(agrunc,2))),pch=NA,lty=1,bty='n',lwd=3)

	dev.off()

	## EPS Version

	# --- Dimensions and Scaling Logic ---
	width_cm <- 8.9
	width_in <- width_cm / 2.54
	sc <- 0.21 # Scaling factor for absolute line widths
	new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)

	file_name <- paste('oco_v10mip_fluxes_trop_vs_et_total_byMod_', starttxt, '_to_', endtxt, '.eps', sep='')

	cairo_ps(filename = file_name,
		 width = width_in,
		 height = width_in,
		 pointsize = new_pointsize)

	par(mar=c(5, 5, 4, 1)+0.1)

	modmeantrop=NULL
	modsdtrop=NULL
	modmeanet=NULL
	modsdet=NULL
	for(mod in mods){
		modmeantrop=c(modmeantrop,mean(tropmeanfluxes$total[tropmeanfluxes$mod==mod]))
		modsdtrop=c(modsdtrop,sd(tropmeanfluxes$total[tropmeanfluxes$mod==mod]))
		modmeanet=c(modmeanet,mean(etmeanfluxes$total[etmeanfluxes$mod==mod]))
		modsdet=c(modsdet,sd(etmeanfluxes$total[etmeanfluxes$mod==mod]))
	}

	plot(etmeanfluxes$total,tropmeanfluxes$total,xlab=expression(paste('Extratropical Total Flux (PgC ',yr^-1,')')),ylab=expression(paste('Tropical Total Flux (PgC ',yr^-1,')')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Tropical vs. Extratropical Total Fluxes'),type='n',axes=F)
	xcorners=c(-1,7,7,-1,-1)
	ycorners=-1*xcorners+agr+c(-1*agrunc,-1*agrunc,agrunc,agrunc,-1*agrunc)
	polygon(xcorners,ycorners,col='gray75',border=NA)
	abline(agr,-1,lwd=3*sc)

	box(lwd=0.5)
	axis(1,cex.axis=1.5,lwd=0.5)
	axis(2,cex.axis=1.5,lwd=0.5)

	modord=order(modsdet,decreasing=T)
	for(mod in mods[modord]){
		expord=order(etmeanfluxes$total[etmeanfluxes$mod==mod])
		exppch=match(etmeanfluxes$exp[etmeanfluxes$mod==mod][expord],exps)+20
		exppch[exppch==26]=21 # added for testing unopt
		points(etmeanfluxes$total[etmeanfluxes$mod==mod][expord],tropmeanfluxes$total[tropmeanfluxes$mod==mod][expord],col=cols[match(mod,mods)],cex=2,lwd=3*sc,pch=exppch,bg=bgcols[match(mod,mods)],type='b') # ,'white')
	}
	legend('topright',modnames,pch=NA,lty=NA,text.col=cols[c(1:length(mods))],bty='n',title='Model:',title.col='black',title.adj=0.5,lwd=1)
	legend('bottomleft',exps,col='black',pch=21:25,lty=NA,lwd=3*sc,pt.cex=1.5,bty='n',pt.bg='white',title='Experiment:',title.col='black',title.adj=0.1)
	legend('top',legend=substitute(AGR == a %+-% b,list(a=round(agr,2),b=round(agrunc,2))),pch=NA,lty=1,bty='n',lwd=3*sc)

	dev.off()


	png(paste('oco_v10mip_fluxes_trop_vs_et_land_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)

	plot(etmeanfluxes$land,tropmeanfluxes$land,xlab='ET Land Flux (PgC/yr)',ylab='Trop Land Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Trop vs. ET Land Fluxes',type='n')
	for(mod in mods){
		points(etmeanfluxes$land[etmeanfluxes$mod==mod],tropmeanfluxes$land[tropmeanfluxes$mod==mod],col=expcols[match(etmeanfluxes$exp[etmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white')
	}

	for(exp in exps){
		expmeantrop=mean(tropmeanfluxes$land[tropmeanfluxes$exp==exp])
		expsdtrop=sd(tropmeanfluxes$land[tropmeanfluxes$exp==exp])
		expmeanet=mean(etmeanfluxes$land[etmeanfluxes$exp==exp])
		expsdet=sd(etmeanfluxes$land[etmeanfluxes$exp==exp])
		points(expmeanet,expmeantrop,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeanet-expsdet,expmeantrop,expmeanet+expsdet,expmeantrop,lwd=4,col=expcols[match(exp,exps)])
		segments(expmeanet,expmeantrop-expsdtrop,expmeanet,expmeantrop+expsdtrop,lwd=4,col=expcols[match(exp,exps)])
	}

	dev.off()


	png(paste('oco_v10mip_fluxes_trop_vs_et_ocean_',starttxt,'_to_',endtxt,'.png',sep=''),width=1200,height=1200,pointsize=24)

	plot(etmeanfluxes$ocean,tropmeanfluxes$ocean,xlab='ET Ocean Flux (PgC/yr)',ylab='Trop Ocean Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Trop vs. ET Ocean Fluxes',type='n')
	for(mod in mods){
		points(etmeanfluxes$ocean[etmeanfluxes$mod==mod],tropmeanfluxes$ocean[tropmeanfluxes$mod==mod],col=expcols[match(etmeanfluxes$exp[etmeanfluxes$mod==mod],exps)],cex=2,lwd=3,pch=match(mod,mods),bg='white')
	}

	for(exp in exps){
		expmeantrop=mean(tropmeanfluxes$ocean[tropmeanfluxes$exp==exp])
		expsdtrop=sd(tropmeanfluxes$ocean[tropmeanfluxes$exp==exp])
		expmeanet=mean(etmeanfluxes$ocean[etmeanfluxes$exp==exp])
		expsdet=sd(etmeanfluxes$ocean[etmeanfluxes$exp==exp])
		points(expmeanet,expmeantrop,cex=2,col=expcols[match(exp,exps)],pch=16)
		segments(expmeanet-expsdet,expmeantrop,expmeanet+expsdet,expmeantrop,lwd=4,col=expcols[match(exp,exps)])
		segments(expmeanet,expmeantrop-expsdtrop,expmeanet,expmeantrop+expsdtrop,lwd=4,col=expcols[match(exp,exps)])
	}

	dev.off()


	pngdim=1500
	pngpts=30

	png('flux_vs_flux_legend_byExp.png',height=pngdim,width=pngdim,pointsize=pngpts)
	plot(c(1:length(mods)),c(1:length(mods)),type='n',xaxt='n',yaxt='n',xlab='',ylab='')
	legend('topleft',mods,col='black',cex=2,pt.lwd=3,pch=c(1:length(mods)),pt.bg='white',bty='n',title='Model:',title.col='black',title.adj=0.1)
	legend('bottomright',exps,col=cols[c(1:length(exps))],cex=2,pt.lwd=3,pch=0,pt.bg='white',bty='n',title='Experiment:',title.col='black',title.adj=0.1)
	dev.off()

	png('flux_vs_flux_legend_byMod.png',height=pngdim,width=pngdim,pointsize=pngpts)
	plot(c(1:length(mods)),c(1:length(mods)),type='n',xaxt='n',yaxt='n',xlab='',ylab='')
	legend('topleft',mods,col=cols[c(1:length(mods))],cex=2,pt.lwd=3,pch=0,pt.bg='white',bty='n',title='Model:',title.col='black',title.adj=0.1)
	legend('bottomright',exps,pch=c(1:length(exps)),col='black',cex=2,pt.lwd=3,pt.bg='white',bty='n',title='Experiment:',title.col='black',title.adj=0.1)
	dev.off()

	# make paper versions
	system(paste('montage oco_v10mip_fluxes_ocean_vs_land_byMod_',starttxt,'_to_',endtxt,'.png oco_v10mip_fluxes_trop_vs_et_total_byExp_',starttxt,'_to_',endtxt,'.png -geometry +1+2 fig1_dipoles_',starttxt,'_to_',endtxt,'.png',sep=''))
	system(paste('montage oco_v10mip_fluxes_ocean_vs_land_byExp_',starttxt,'_to_',endtxt,'.png oco_v10mip_fluxes_trop_vs_et_total_byMod_',starttxt,'_to_',endtxt,'.png -geometry +1+2 figS1_dipoles_',starttxt,'_to_',endtxt,'.png',sep=''))
	system(paste('cp oco_v10mip_fluxes_ocean_vs_land_byMod_',starttxt,'_to_',endtxt,'.eps fig1a.eps',sep=''))
	system(paste('cp oco_v10mip_fluxes_trop_vs_et_total_byExp_',starttxt,'_to_',endtxt,'.eps fig1b.eps',sep=''))
	system(paste('cp oco_v10mip_fluxes_ocean_vs_land_byExp_',starttxt,'_to_',endtxt,'.eps figS1a.eps',sep=''))
	system(paste('cp oco_v10mip_fluxes_trop_vs_et_total_byMod_',starttxt,'_to_',endtxt,'.eps figS1b.eps',sep=''))

} # end of function
