# program to read in OCO2 MIP model flux files and aggregate and plot

### move fossfile to common LAC data dir and give original source
### still need separate OU dir?

library('ncdf4')
library('RColorBrewer')
library('abind')
cols=brewer.pal(5,'Set1')[c(1,3,2,4,5)]

## Options:
exps=c('IS','OG','LNLG','LNLGOGIS','LNLGIS','unopt') # no LoFI Prior
mods=c('AMES','BAKER','CAMS','CMS-FLUX','COLA','CSU','CT','NIES','OU','TM5-4DVAR','UT','WEIR','WOMBAT') # all 14 except JHU
# JHU ocean fluxes are all identical and ~ 1 PgC/yr high

fossdir='../LOCALDATA'
fossfile='SFCO2_FF.OCO2-MIP.v2020.1.1x1.19860101-20201231.nc'
inclmthetae=F

lonlim=c(-180,180) # for flux calculations

cacalclist=c( # slat, nlat, limtxt, plot
20,90,'20N-90N',T,
-20,20,'20S-20N',T,
-90,-20,'90S-20S',T,
-90,-80,'90S-80S',F,
-80,-70,'80S-70S',F,
-70,-60,'70S-60S',F,
-60,-50,'60S-50S',F,
-50,-40,'50S-40S',F,
-40,-30,'40S-30S',F,
-30,-20,'30S-20S',F,
-20,-10,'20S-10S',F,
-10,0,'10S-EQ',F,
0,10,'EQ-10N',F,
10,20,'10N-20N',F,
20,30,'20N-30N',F,
30,40,'30N-40N',F,
40,50,'40N-50N',F,
50,60,'50N-60N',F,
60,70,'60N-70N',F,
70,80,'70N-80N',F,
80,90,'80N-90N',F,
20,90,'20X-90X',T # lump together SET and NET
)

caspecs=matrix(cacalclist,byrow=T,ncol=4)

# Process fluxes:

# open fossil fuel flux file (OCO v10 MIP product from APO fwd repo)
fossnc=nc_open(paste(fossdir,'/',fossfile,sep=''))
fossfluxdaily=ncvar_get(fossnc,'SFCO2_FF') # float SFCO2_FF[lon,lat,time], units: mol/m^2/s
fossarea=ncvar_get(fossnc,'area') # float area[lon,lat], units: m^2
# lat and lon grid same as ocean and land fluxes below
# daily resolution so have to aggregate to monthly
# also need to convert to PgC
# also need to trim to 2015:2020
# assign a month index corresponding to land/ocean fluxes
tmcomp=ncvar_get(fossnc,'time_components') # float time_components[n_time_components,time], long_name: time components (year, month, day, hour, min, sec)
fossyear=tmcomp[1,]
fossmondaily=tmcomp[2,]
fossmon=rep(0,dim(fossfluxdaily)[3])
fossmon[fossyear>=2015&fossyear<=2020]=fossmondaily[fossyear>=2015&fossyear<=2020]+(fossyear[fossyear>=2015&fossyear<=2020]-2015)*12 # 1-72 for Jan 2015 - Dec 2020


print('looping on experiments in:')
print(exps)
## subdirs are allcaps
#> unique(opncsubs)
# [1] "CT"        "CAMS"      "UT"        "Baker"     "WOMBAT"    "CSU"
# [7] "TM5-4DVAR" "OU"        "Ames"      "COLA"      "CMS-Flux"  "LoFI"
#[13] "M2CC"      "JHU"       "CSU-NEE"   "NIES"

print('looping on models in:')
print(mods)

# loop on experiment
for(exp in exps){

        print(exp)
	expname=exp
	if(expname=='unopt') expname='Prior'

        # loop on models
        for(mod in mods){

		print(mod)

                # skip WEIR for all but IS
                if(mod!='WEIR'|exp=='IS'){

			if(!file.exists(mod)) system(paste('mkdir',mod))
			if(!file.exists(paste(mod,'/',exp,sep=''))) system(paste('mkdir ',mod,'/',exp,sep=''))

                        fluxname=mod
			# V10: No ATom file for CSU, UT, JHU (also missing 3 NIES)
			fluxdir='../LOCALDATA/V10MIP'
			if(fluxname=='CMS-FLUX') fluxname='CMS-Flux'; if(fluxname=='AMES') fluxname='Ames'; if(fluxname=='BAKER') fluxname='Baker';
			if(fluxname=='WEIR') fluxname='LoFI'
			### clean up after testing
#			if(fluxname!='OU'){
				flux=nc_open(paste(fluxdir,'/',expname,'/',fluxname,'_gridded_fluxes_',expname,'.nc4',sep=''))
#			} else {
#				# new fixed files
#				flux=nc_open(paste(fluxdir,'/fixedOUfluxes/',fluxname,'_gridded_fluxes_',expname,'.nc4',sep=''))
#			}
			fmon=rep(seq(1,12),6) 
			year=rep(2015:2020,each=12)
			moddectime=year+(fmon-0.5)/12
                        flat=ncvar_get(flux,'latitude') # box centers
                        flon=ncvar_get(flux,'longitude') # box centers
                        landflux=ncvar_get(flux,'land') # 360 x 180 x 48 (v9) or 72 (v10)
                        oceanflux=ncvar_get(flux,'ocean')
			fossflux=array(NA,dim=dim(oceanflux))
			### clean up after testing
#			if(fluxname!='OU'){
			       netflux=ncvar_get(flux,'net')
#			} else if(fluxname=='OU'){
#			       netflux=landflux+oceanflux
#			}
                        # need to convert from "gC per m2 per year" to PgC/yr
			# area of sphere = 4*pi*rearth^2, area of zone = 2*pi*rearth*h where h from equator = rearth*sin(lat)
			rearth=6.371E6 # m
			southedge=flat-0.5
			northedge=flat+0.5
			surfarea=(2*pi*rearth*(rearth*sin(northedge*pi/180))-2*pi*rearth*(rearth*sin(southedge*pi/180)))/360
			for(mon in c(1:72)){
				landflux[,,mon]=t(t(landflux[,,mon])*as.vector(surfarea))/1E15
				oceanflux[,,mon]=t(t(oceanflux[,,mon])*as.vector(surfarea))/1E15
				netflux[,,mon]=t(t(netflux[,,mon])*as.vector(surfarea))/1E15
				fossflux[,,mon]=t(t(apply(fossfluxdaily[,,fossmon==mon],c(1,2),mean))*as.vector(surfarea))*12/1E15*86400*365 # converts to PgC/yr averaged over month
			}

			landfluxsave=landflux
			oceanfluxsave=oceanflux
			netfluxsave=netflux
			fossfluxsave=fossflux
			
			for(i in c(1:nrow(caspecs))){

				landflux=landfluxsave
				oceanflux=oceanfluxsave
				netflux=netfluxsave
				fossflux=fossfluxsave

				print(caspecs[i,])
				latlim=as.numeric(caspecs[i,1:2])
				limtxt=caspecs[i,3]

				# mask for region
				if(!grepl('MThetae',limtxt)){
					
					if(grepl('X',limtxt)){ # lump together SH and NH
						landflux[flon<lonlim[1]|flon>=lonlim[2],,]=0; landflux[,abs(flat)<latlim[1]|abs(flat)>=latlim[2],]=0
						oceanflux[flon<lonlim[1]|flon>=lonlim[2],,]=0; oceanflux[,abs(flat)<latlim[1]|abs(flat)>=latlim[2],]=0
					} else {
						landflux[flon<lonlim[1]|flon>=lonlim[2],,]=0; landflux[,flat<latlim[1]|flat>=latlim[2],]=0
						oceanflux[flon<lonlim[1]|flon>=lonlim[2],,]=0; oceanflux[,flat<latlim[1]|flat>=latlim[2],]=0
					}
					# sum region
					reglandflux=apply(landflux,3,sum)
					regoceanflux=apply(oceanflux,3,sum)
					if(grepl('X',limtxt)){ # lump together SH and NH
						netflux[flon<lonlim[1]|flon>=lonlim[2],,]=0; netflux[,abs(flat)<latlim[1]|abs(flat)>=latlim[2],]=0
						fossflux[flon<lonlim[1]|flon>=lonlim[2],,]=0; fossflux[,abs(flat)<latlim[1]|abs(flat)>=latlim[2],]=0
					} else {
						netflux[flon<lonlim[1]|flon>=lonlim[2],,]=0; netflux[,flat<latlim[1]|flat>=latlim[2],]=0
						fossflux[flon<lonlim[1]|flon>=lonlim[2],,]=0; fossflux[,flat<latlim[1]|flat>=latlim[2],]=0
					}
					regnetflux=apply(netflux,3,sum)
					regfossflux=apply(fossflux,3,sum)

				} else { # use MThetae
					if(grepl('NH',limtxt)){ 
						mtemask=(mtenh/1E16>latlim[1]&mtenh/1E16<latlim[2])*1
						mtemask=abind(mtemask*0,mtemask,along=2) # set SH to 0
					} else if(grepl('SH',limtxt)){ 
						mtemask=(mtesh/1E16>latlim[1]&mtesh/1E16<latlim[2])*1
						mtemask=abind(mtemask,mtemask*0,along=2) # set NH to 0
					}
					# loop through by day adding up
					reglandflux=rep(0,dim(landflux)[3])
					regoceanflux=rep(0,dim(oceanflux)[3])
					regnetflux=rep(0,dim(netflux)[3])
					regfossflux=rep(0,dim(fossflux)[3])
					allmon=(mtedt$year+1900-2015)*12+mtedt$mon+1
                                        for(j in c(1:length(mtedt))){
                                                mon=allmon[j]
                                                reglandflux[mon]=reglandflux[mon]+sum(landflux[,,mon]*mtemask[,,j]/ydays[j]) # fluxesin PgC/yr, adding PgC/day
                                                regoceanflux[mon]=regoceanflux[mon]+sum(oceanflux[,,mon]*mtemask[,,j]/ydays[j])
                                                regnetflux[mon]=regnetflux[mon]+sum(netflux[,,mon]*mtemask[,,j]/ydays[j])
                                                regfossflux[mon]=regfossflux[mon]+sum(fossflux[,,mon]*mtemask[,,j]/ydays[j])
                                        }
					reglandflux=reglandflux*12 # back to PgC/yr
					regoceanflux=regoceanflux*12 
					regnetflux=regnetflux*12 
					regfossflux=regfossflux*12

				}

				# write out flux timeseries
				write('year month land ocean net foss',paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''))
				write(rbind(year,fmon,reglandflux,regoceanflux,regnetflux,regfossflux),ncol=6,paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''),append=T)

				# aggregate by month
				climlandflux=aggregate(reglandflux,by=list(mon=fmon),mean) # X-yr mean
				climoceanflux=aggregate(regoceanflux,by=list(mon=fmon),mean) # X-yr mean
				climnetflux=aggregate(regnetflux,by=list(mon=fmon),mean) # X-yr mean
				climfossflux=aggregate(regfossflux,by=list(mon=fmon),mean) # X-yr mean

				if(caspecs[i,4]){

					## plot fluxes
					print('plotting')

					png(paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.png',sep=''),height=1200,width=1800,pointsize=30)
					par(mar=c(5,5,4,2)+0.1)

					ylm=range(c(reglandflux,regoceanflux,regnetflux,regfossflux))
					plot(moddectime,reglandflux,type='b',xlab='Year',ylab='PgC/yr',main=paste(mod,' ',exp,'  flux (',latlim[1],' to ',latlim[2],', ',lonlim[1],' to ',lonlim[2],')',sep=''),cex.main=1.2,cex.lab=1.2,pch=1,col=cols[2],cex=1.5,lwd=2,ylim=ylm)
					points(moddectime,regoceanflux,type='b',pch=1,col=cols[3],cex=1.5,lwd=2)
					points(moddectime,regnetflux,type='b',pch=1,col=cols[4],cex=1.5,lwd=2)
					points(moddectime,regfossflux,type='b',pch=1,col=cols[5],cex=1.5,lwd=2)
					legend('bottomright',c('Land','Ocean','Net','Foss'),col=c(cols[2],cols[3],cols[4],cols[5]),pch=1,pt.lwd=2,cex=1.5,lty=1)
					dev.off()

					png(paste(mod,'/',exp,'/flux_',limtxt,'_seascycle.png',sep=''),height=1200,width=1800,pointsize=30)
					par(mar=c(5,5,4,2)+0.1)
					
					ylm=range(c(climlandflux$x,climoceanflux$x,climnetflux$x,climfossflux$x))
					plot(seq(0.5,11.5),climlandflux$x,type='n',xlab='Year',ylab='PgC',main=paste(mod,' ',exp,' flux (',latlim[1],' to ',latlim[2],', ',lonlim[1],' to ',lonlim[2],')',sep=''),cex.main=1.2,cex.lab=1.2,axes=F,ylim=ylm)
					box()
					axis(2,cex.axis=1.2)
					axis(1,at=c(0:12),labels=F,cex.axis=1.5)
					axis(1,seq(0.5,11.5),labels=c('J','F','M','A','M','J','J','A','S','O','N','D'),cex.axis=1.3,tick=F)
					points(seq(0.5,11.5),climlandflux$x,pch=1,col=cols[2],cex=1.5,lwd=2,type='b')
					points(seq(0.5,11.5),climoceanflux$x,pch=1,col=cols[3],cex=1.5,lwd=2,type='b')
					points(seq(0.5,11.5),climnetflux$x,pch=1,col=cols[4],cex=1.5,lwd=2,type='b')
					points(seq(0.5,11.5),climfossflux$x,pch=1,col=cols[5],cex=1.5,lwd=2,type='b')
					legend('bottomright',c('Land','Ocean','Net','Foss'),col=c(cols[2],cols[3],cols[4],cols[5]),pch=1,pt.lwd=2,cex=1.5,lty=1)
					dev.off()

				}

			} # loop on caspecs

                } # skip some

        } # loop on model

} # loop on experiment

## now loop on caspecs, plotting flux timeseries
dummyflux=read.table(paste('AMES/IS/flux_20N-90N_timeseries.txt',sep=''),header=T)
dummyflux$total=dummyflux$land+dummyflux$ocean+dummyflux$foss
dummyfluxdt=as.POSIXlt(ISOdate(dummyflux$year,dummyflux$month,15,0),tz='UTC') # midnight on 15th of each month

ylm=c(-35,25)
cols=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2')) # 17 diff colors
for(i in c(1:nrow(caspecs))){

	if(caspecs[i,4]){

		limtxt=caspecs[i,3]
		latlim=as.numeric(caspecs[i,1:2])
		png(paste('flux_ts_',limtxt,'_byexp.png',sep=''),height=1200,width=1800,pointsize=30)
		par(mfrow=c(3,2))

		# loop on experiment
		for(exp in exps){

			expmean=rep(0,nrow(dummyflux))
			nmean=0
			plot(dummyfluxdt,dummyflux$total,type='n',xlab='Year',ylab='PgC/yr',main=paste(exp,' Fluxes (',latlim[1],' to ',latlim[2],', ',lonlim[1],' to ',lonlim[2],')',sep=''),cex.main=1.2,cex.lab=1.2,pch=1,col=cols[2],cex=1.5,lwd=2,ylim=ylm)

			# loop on models
			for(mod in mods){

				if(mod!='WEIR'|exp=='IS'){

					monflux=read.table(paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''),header=T)
					monflux$total=monflux$land+monflux$ocean+monflux$foss
					monfluxdt=as.POSIXlt(ISOdate(monflux$year,monflux$month,15,0),tz='UTC') # midnight on 15th of each month
					points(monfluxdt,monflux$total,type='b',pch=1,col=cols[match(mod,mods)],cex=1.5,lwd=2)
					expmean=expmean+monflux$total
					nmean=nmean+1

				} # skip some

			} # loop on model
			assign(paste(exp,'mean',sep=''),expmean/nmean)

		} # loop on experiment

		dev.off()

		png(paste('flux_ts_',limtxt,'_expmean.png',sep=''),height=1200,width=1800,pointsize=30)

		plot(dummyfluxdt,dummyflux$total,type='n',xlab='Year',ylab='PgC/yr',main=paste('Experiment Mean Fluxes (',latlim[1],' to ',latlim[2],', ',lonlim[1],' to ',lonlim[2],')',sep=''),cex.main=1.2,cex.lab=1.2,pch=1,col=cols[2],cex=1.5,lwd=2,ylim=ylm)
		for(exp in exps){
			points(dummyfluxdt,get(paste(exp,'mean',sep='')),type='b',pch=1,col=cols[match(exp,exps)],cex=1.5,lwd=2)
		}

		dev.off()
	
	}

} # loop on caspecs

png('flux_ts_legend_byMod.png',height=1200,width=1800,pointsize=30)
plot(dummyfluxdt,dummyflux$total,type='n',xaxt='n',yaxt='n',xlab='',ylab='')
legend('topleft',mods,col=cols[c(1:length(mods))],cex=2,pt.lwd=3,pch=0,pt.bg='white',bty='n',title='Model:',title.col='black',title.adj=0.1)
dev.off()
