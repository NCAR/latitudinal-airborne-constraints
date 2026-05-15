# program to read in flux_proc.r output and make all-experiment plot vs latitude


library('ncdf4')
library('RColorBrewer')
expcols=c(brewer.pal(5,'Set1'),'black')

## Options:
exps=c('IS','OG','LNLG','LNLGOGIS','LNLGIS') # ,'unopt')
mods=c('AMES','BAKER','CAMS','CMS-FLUX','COLA','CSU','CT','NIES','OU','TM5-4DVAR','UT','WEIR','WOMBAT') # all 14 except JHU

cacalclist=c( # slat, nlat, limtxt
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
caspecs=matrix(cacalclist,byrow=T,ncol=3)

# loop on caspecs, calculating expmeans vs lat for 6/2016-5/2018
start=ISOdate(2016,6,1,tz='UTC') ## 24 months ending the same month as ATom, with 2-months lead before ATom
end=ISOdate(2018,5,31,tz='UTC')
lats=seq(-85,85,10)

# loop on experiment
globalmeanfluxes=NULL
globalmeanfluxesfulltime=NULL
reglim=c(20,90)
allmeanbylat=NULL
allmeanbylatfulltime=NULL
allmeanbylatland=NULL
allmeanbylatocean=NULL

for(exp in exps){

	expmean=NULL
	expmeanocean=NULL
	expmeanland=NULL

	# loop on models
	for(mod in mods){

		# skip WEIR for all but IS
                if(mod!='WEIR'|exp=='IS'){ # &exp!='unopt')){

			modmean=NULL
			modmeanfulltime=NULL
			modmeanocean=NULL
			modmeanland=NULL
			globalmean=0
			globalmeanfulltime=0

			for(i in c(1:nrow(caspecs))){
				limtxt=caspecs[i,3]
				latlim=as.numeric(caspecs[i,1:2]) ### get rid of latlim?

				if(file.exists(paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''))){

					monflux=read.table(paste(mod,'/',exp,'/flux_',limtxt,'_timeseries.txt',sep=''),header=T)
					monflux$total=monflux$land+monflux$ocean+monflux$foss
					monfluxdt=as.POSIXlt(ISOdate(monflux$year,monflux$month,15,0),tz='UTC') # midnight on 15th of each month
					modmean=c(modmean,mean(monflux$total[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0]))
					modmeanfulltime=c(modmeanfulltime,mean(monflux$total))
					globalmean=globalmean+mean(monflux$total[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0]) # adding up over all latbins
					globalmeanfulltime=globalmeanfulltime+mean(monflux$total) # adding up over all latbins
					modmeanocean=c(modmeanocean,mean(monflux$ocean[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0]))
					modmeanland=c(modmeanland,mean(monflux$land[difftime(monfluxdt,start)>0&difftime(monfluxdt,end)<0]))

				}

			} # loop on caspecs

			allmeanbylat=cbind(allmeanbylat,c(exp,mod,round(modmean,3)))
			allmeanbylatfulltime=cbind(allmeanbylatfulltime,c(exp,mod,round(modmeanfulltime,3)))
			allmeanbylatland=cbind(allmeanbylatland,c(exp,mod,round(modmeanland,3)))
			allmeanbylatocean=cbind(allmeanbylatocean,c(exp,mod,round(modmeanocean,3)))
			globalmeanfluxes=rbind(globalmeanfluxes,c(exp,mod,round(globalmean,3)))
			globalmeanfluxesfulltime=rbind(globalmeanfluxesfulltime,c(exp,mod,round(globalmeanfulltime,3))) ### for Bianca
			expmean=rbind(expmean,modmean)
			expmeanocean=rbind(expmeanocean,modmeanocean)
			expmeanland=rbind(expmeanland,modmeanland)

		} # if not WEIR or IS

	} # loop on model
	expmean=apply(expmean,2,mean)
	assign(paste(exp,'mean',sep=''),expmean)
	expmeanocean=apply(expmeanocean,2,mean)
	assign(paste(exp,'meanocean',sep=''),expmeanocean)
	expmeanland=apply(expmeanland,2,mean)
	assign(paste(exp,'meanland',sep=''),expmeanland)

} # loop on experiment

write('exp mod flux','oco2_v10mip_global_total_fluxes_2015-2020.txt')
write(t(globalmeanfluxesfulltime),'oco2_v10mip_global_total_fluxes_2015-2020.txt',ncol=3,append=T)

write('exp mod flux','oco2_v10mip_global_total_fluxes_201606-201805.txt')
write(t(globalmeanfluxes),'oco2_v10mip_global_total_fluxes_201606-201805.txt',ncol=3,append=T)

write(c('exp','mod',lats),'oco2_v10mip_total_fluxes_bylat_2015-2020.txt',ncol=nrow(allmeanbylatfulltime))
write(allmeanbylat,'oco2_v10mip_total_fluxes_bylat_2015-2020.txt',ncol=nrow(allmeanbylatfulltime),append=T)

write(c('exp','mod',lats),'oco2_v10mip_total_fluxes_bylat_201606-201805.txt',ncol=nrow(allmeanbylat))
write(allmeanbylat,'oco2_v10mip_total_fluxes_bylat_201606-201805.txt',ncol=nrow(allmeanbylat),append=T)

write(c('exp','mod',lats),'oco2_v10mip_land_fluxes_bylat_201606-201805.txt',ncol=nrow(allmeanbylatland))
write(allmeanbylatland,'oco2_v10mip_land_fluxes_bylat_201606-201805.txt',ncol=nrow(allmeanbylatland),append=T)

write(c('exp','mod',lats),'oco2_v10mip_ocean_fluxes_bylat_201606-201805.txt',ncol=nrow(allmeanbylatocean))
write(allmeanbylatocean,'oco2_v10mip_ocean_fluxes_bylat_201606-201805.txt',ncol=nrow(allmeanbylatocean),append=T)


png('oco_v10mip_fluxes_ann_exp_means_total.png',width=1800,height=1200,pointsize=24)
ylm=range(c(ISmean,LNLGmean,LNLGISmean,OGmean,LNLGOGISmean))
plot(lats,ISmean,type='n',xlab='Latitude (deg N)',ylab='Total Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP 10-degree Total Fluxes',ylim=ylm)
abline(h=0)
for(exp in exps){
        points(lats,get(paste(exp,'mean',sep='')),type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
}
legend('topleft',exps,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()


png('oco_v10mip_fluxes_ann_exp_means_intfrompole_total.png',width=1800,height=1200,pointsize=24)
for(exp in exps){
	tmp=get(paste(exp,'mean',sep=''))
	tmp[lats<0]=cumsum(tmp[lats<0])
	tmp[lats>0]=rev(cumsum(rev(tmp[lats>0])))
	assign(paste(exp,'int',sep=''),tmp)
}
ylm=range(c(ISint,LNLGint,LNLGISint,OGint,LNLGOGISint))
plot(lats,ISint,type='n',xlab='Latitude (deg N)',ylab='Total Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Total Fluxes Integrated from Poles',ylim=ylm)
abline(h=0)
for(exp in exps){
        points(lats[lats<0],get(paste(exp,'int',sep=''))[lats<0],type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
        points(lats[lats>0],get(paste(exp,'int',sep=''))[lats>0],type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
}
legend('topleft',exps,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()


png('oco_v10mip_fluxes_ann_exp_means_ocean.png',width=1800,height=1200,pointsize=24)
ylm=range(c(ISmeanocean,LNLGmeanocean,LNLGISmeanocean,OGmeanocean,LNLGOGISmeanocean))
plot(lats,ISmeanocean,type='n',xlab='Latitude (deg N)',ylab='Ocean Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP 10-degree Ocean Fluxes',ylim=ylm,xaxp=c(-90,90,6))
abline(h=0)
for(exp in exps){
        points(lats,get(paste(exp,'meanocean',sep='')),type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
}
legend('topleft',exps,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()


png('oco_v10mip_fluxes_ann_exp_means_intfrompole_ocean.png',width=1800,height=1200,pointsize=24)
for(exp in exps){
        tmp=get(paste(exp,'meanocean',sep=''))
        tmp[lats<0]=cumsum(tmp[lats<0])
        tmp[lats>0]=rev(cumsum(rev(tmp[lats>0])))
        assign(paste(exp,'intocean',sep=''),tmp)
}
ylm=range(c(ISintocean,LNLGintocean,LNLGISintocean,OGintocean,LNLGOGISintocean))
plot(lats,ISintocean,type='n',xlab='Latitude (deg N)',ylab='Ocean Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Ocean Fluxes Integrated from Poles',ylim=ylm,xaxp=c(-90,90,6))
abline(h=0)
for(exp in exps){
        points(lats[lats<0],get(paste(exp,'intocean',sep=''))[lats<0],type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
        points(lats[lats>0],get(paste(exp,'intocean',sep=''))[lats>0],type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
}
legend('bottomleft',exps,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()


png('oco_v10mip_fluxes_ann_exp_means_land.png',width=1800,height=1200,pointsize=24)
ylm=range(c(ISmeanland,LNLGmeanland,LNLGISmeanland,OGmeanland,LNLGOGISmeanland))
plot(lats,ISmeanland,type='n',xlab='Latitude (deg N)',ylab='Land Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP 10-degree Land Fluxes',ylim=ylm,xaxp=c(-90,90,6))
abline(h=0)
for(exp in exps){
        points(lats,get(paste(exp,'meanland',sep='')),type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
}
legend('topleft',exps,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()


png('oco_v10mip_fluxes_ann_exp_means_intfrompole_land.png',width=1800,height=1200,pointsize=24)
for(exp in exps){
        tmp=get(paste(exp,'meanland',sep=''))
        tmp[lats<0]=cumsum(tmp[lats<0])
        tmp[lats>0]=rev(cumsum(rev(tmp[lats>0])))
        assign(paste(exp,'intland',sep=''),tmp)
}
ylm=range(c(ISintland,LNLGintland,LNLGISintland,OGintland,LNLGOGISintland))
plot(lats,ISintland,type='n',xlab='Latitude (deg N)',ylab='Land Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Land Fluxes Integrated from Poles',ylim=ylm,xaxp=c(-90,90,6))
abline(h=0)
for(exp in exps){
        points(lats[lats<0],get(paste(exp,'intland',sep=''))[lats<0],type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
        points(lats[lats>0],get(paste(exp,'intland',sep=''))[lats>0],type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
}
legend('bottomleft',exps,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()


png('oco_v10mip_fluxes_ann_exp_means_ocean_minus_land.png',width=1800,height=1200,pointsize=24)
ylm=range(c(ISmeanocean-ISmeanland,LNLGmeanocean-LNLGmeanland,LNLGISmeanocean-LNLGISmeanland,OGmeanocean-OGmeanland,LNLGOGISmeanocean-LNLGOGISmeanland))
plot(lats,ISmeanocean-ISmeanland,type='n',xlab='Latitude (deg N)',ylab='Ocean - Land Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP 10-degree Ocean minus Land Fluxes',ylim=ylm,xaxp=c(-90,90,6))
abline(h=0)
for(exp in exps){
        points(lats,get(paste(exp,'meanocean',sep=''))-get(paste(exp,'meanland',sep='')),type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
}
legend('topleft',exps,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()


png('oco_v10mip_fluxes_ann_exp_means_intfrompole_ocean_minus_land.png',width=1800,height=1200,pointsize=24)
ylm=range(c(ISintocean-ISintland,LNLGintocean-LNLGintland,LNLGISintocean-LNLGISintland,OGintocean-OGintland,LNLGOGISintocean-LNLGOGISintland))
plot(lats,ISintocean-ISintland,type='n',xlab='Latitude (deg N)',ylab='Ocean - Land Flux (PgC/yr)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='OCO-2 V10MIP Ocean minus Land Fluxes Integrated from Poles',ylim=ylm,xaxp=c(-90,90,6))
abline(h=0)
for(exp in exps){
        points(lats[lats<0],get(paste(exp,'intocean',sep=''))[lats<0]-get(paste(exp,'intland',sep=''))[lats<0],type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
        points(lats[lats>0],get(paste(exp,'intocean',sep=''))[lats>0]-get(paste(exp,'intland',sep=''))[lats>0],type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(exps==exp)])
}
legend('topleft',exps,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()
