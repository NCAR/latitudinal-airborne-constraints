# program to plot column-mean latitudinal concentration gradients

modeldir='V10MIP'
obsdir='OBS/BINFIT'


colval=read.table(paste(obsdir,'/Col_Values_CO2_OP_10x100.txt',sep=''),header=T)
colnames(colval)[4]='co2'
if(is.element('slc',colnames(colval))){
	aggby='campslc'
} else if(is.element('flt',colnames(colval))){
	aggby='campflt'
} else if(is.element('yday',colnames(colval))){
	aggby='day'
}
anncol=read.table(paste(obsdir,'/Annual_Mean_CO2_OP_10x100_Column.txt',sep=''))
colnames(anncol)=c('lat','co2')
obsanncol=anncol
annsurf=read.table(paste(obsdir,'/Annual_Mean_CO2_OP_10x100_Surface.txt',sep=''))
colnames(annsurf)=c('lat','co2')
obsannsurf=annsurf

gradco=900

annmat=read.table(paste(obsdir,'/Annual_Mean_CO2_OP_10x100.txt',sep=''),skip=2)
# lats (StoN) as rows, prs (350->950) as columns
latvals=scan(paste(obsdir,'/Annual_Mean_CO2_OP_10x100.txt',sep=''),nlines=1)
prsvals=scan(paste(obsdir,'/Annual_Mean_CO2_OP_10x100.txt',sep=''),skip=1,nlines=1)
if(sum(prsvals>gradco)>1){
	anngrad=cbind(latvals,apply(annmat[,prsvals>gradco],1,mean)-apply(annmat[,prsvals<gradco],1,mean))
} else {
	anngrad=cbind(latvals,annmat[,prsvals>gradco]-apply(annmat[,prsvals<gradco],1,mean))
}
anngrad=data.frame(anngrad)
colnames(anngrad)=c('lat','co2')
obsanngrad=anngrad


library('RColorBrewer')
library('colorspace')
cols=brewer.pal(4,'Set1')
#modcols=brewer.pal(12,'Paired'); modcols[11]='Gold'
modcols=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2')) # 17 diff colors
expcols=brewer.pal(5,'Set1')
pchs=c(6,2) # (0 = Southbound, 1 = Northbound) 

png(paste(modeldir,'/ATom_partial_columns.png',sep=''),width=1800,height=1200,pointsize=24)

plot(colval$lat,colval$co2,type='n',xlab='Latitude (deg N)',ylab='Delta CO2 (- MLO trend)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='ATom 1000-300 mb 10-degree Partial Column Means')

for(camp in c(21:24)){

	if(aggby=='campslc'){
		for(slc in c(0,1)){
			points(colval$lat[colval$camp==camp&colval$slc==slc],colval$co2[colval$camp==camp&colval$slc==slc],type='b',col=cols[camp-20],pch=pchs[slc+1],cex=1,lwd=4)
		}
	} else {
		points(colval$lat[colval$camp==camp],colval$co2[colval$camp==camp],type='b',col=cols[camp-20],pch=16,cex=1,lwd=4)

	}

}
points(anncol$lat[is.element(anncol$lat,colval$lat)],anncol$co2[is.element(anncol$lat,colval$lat)],type='b',pch=16,cex=1.5,lwd=6,col='black')

if(aggby=='campslc'){
	legend('topleft',c('ATom 1 - Aug. 2016','Atom 2 - Feb. 2017','ATom 3 - Oct. 2017','ATom 4 - May 2018','Pacific - Southbound','Atlantic - Northbound','Annual Mean'),col=c(cols,'gray35','gray35','black'),pch=c(rep(11,4),6,2,16),pt.cex=c(rep(1,6),1.5),pt.lwd=2,lwd=c(rep(4,4),NA,NA,6),cex=1.5)
} else {
	legend('topleft',c('ATom 1 - Aug. 2016','Atom 2 - Feb. 2017','ATom 3 - Oct. 2017','ATom 4 - May 2018','Annual Mean'),col=c(cols,'black'),pch=c(rep(16,4),16),pt.cex=c(rep(1,4),1.5),pt.lwd=2,lwd=c(rep(4,4),6),cex=1.5)
}

dev.off()


# V10MIP ATom ObsPack data already merged in with flight data by:
# V10MIP/model_aircraft_merge.r
# which writes files ATom_Merged_Model_Output.txt in subdirectories

## Options:
models=c('AMES','BAKER','CAMS','CMS-FLUX','COLA','CT','OU','TM5-4DVAR','UT','WOMBAT') ### add CSU
experiments=c('IS','LNLG','LNLGIS','LNLGOGIS','OG') # ordered to minimize lines doubling back
expcols=expcols[c(1,3,5,4,2)]

for(exp in experiments){
        print(exp)

	expsetcol=NULL
	expsetgrad=NULL
	expsetsurf=NULL
        for(mod in models){

		if(((exp!='IS'&exp!='LNLG')|mod!='NIES')&(mod!='WEIR'|exp=='IS')){ # no NIES IS or LNLG, and only IS for WEIR

			print(mod)

			colval=read.table(paste(modeldir,'/',mod,'/',exp,'/Col_Values_CO2_OP_10x100.txt',sep=''),header=T) 
			colnames(colval)[4]='co2'
			anncol=read.table(paste(modeldir,'/',mod,'/',exp,'/Annual_Mean_CO2_OP_10x100_Column.txt',sep='')) 
			colnames(anncol)=c('lat','co2')
			annsurf=read.table(paste(modeldir,'/',mod,'/',exp,'/Annual_Mean_CO2_OP_10x100_Surface.txt',sep='')) 
			colnames(annsurf)=c('lat','co2')
			annmat=read.table(paste(modeldir,'/',mod,'/',exp,'/Annual_Mean_CO2_OP_10x100.txt',sep=''),skip=2) # from parent directory binning by flt
			latvals=scan(paste(modeldir,'/',mod,'/',exp,'/Annual_Mean_CO2_OP_10x100.txt',sep=''),nlines=1)
			prsvals=scan(paste(modeldir,'/',mod,'/',exp,'/Annual_Mean_CO2_OP_10x100.txt',sep=''),skip=1,nlines=1)
			if(sum(prsvals>gradco)>1){
				anngrad=cbind(latvals,apply(annmat[,prsvals>gradco],1,mean)-apply(annmat[,prsvals<gradco],1,mean))
			} else {
				anngrad=cbind(latvals,annmat[,prsvals>gradco]-apply(annmat[,prsvals<gradco],1,mean))
			}
			anngrad=data.frame(anngrad)
			colnames(anngrad)=c('lat','co2')
			if(is.null(expsetcol)){
				expsetcol=anncol
				expsetgrad=anngrad
				expsetsurf=annsurf
			} else {
				expsetcol=cbind(expsetcol,anncol$co2)
				expsetgrad=cbind(expsetgrad,anngrad$co2)
				expsetsurf=cbind(expsetsurf,annsurf$co2)
			}
			colnames(expsetcol)[ncol(expsetcol)]=mod
			colnames(expsetgrad)[ncol(expsetgrad)]=mod
			colnames(expsetsurf)[ncol(expsetsurf)]=mod
	
			png(paste(modeldir,'/',mod,'/',exp,'/',mod,'_',exp,'_ATom_partial_columns.png',sep=''),width=1800,height=1200,pointsize=24)
			plot(colval$lat,colval$co2,type='n',xlab='Latitude (deg N)',ylab='Delta CO2 (- MLO trend)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=paste(mod,' ',exp,' ATom 1000-300 mb 10-degree Partial Column Means',sep=''))
			for(camp in c(21:24)){
				if(aggby=='campslc'){
			                for(slc in c(0,1)){
			                	points(colval$lat[colval$camp==camp&colval$slc==slc],colval$co2[colval$camp==camp&colval$slc==slc],type='b',col=cols[camp-20],pch=pchs[slc+1],cex=1,lwd=4)
                			}
			        } else {
			                points(colval$lat[colval$camp==camp],colval$co2[colval$camp==camp],type='b',col=cols[camp-20],pch=16,cex=1,lwd=4)
			        }
			}
			points(anncol$lat[is.element(anncol$lat,colval$lat)],anncol$co2[is.element(anncol$lat,colval$lat)],type='b',pch=16,cex=1.5,lwd=6,col='black')
			if(aggby=='campslc'){
				legend('topleft',c('ATom 1 - Aug. 2016','Atom 2 - Feb. 2017','ATom 3 - Oct. 2017','ATom 4 - May 2018','Pacific - Southbound','Atlantic - Northbound','Annual Mean'),col=c(cols,'gray35','gray35','black'),pch=c(rep(11,4),6,2,16),pt.cex=c(rep(1,6),1.5),pt.lwd=2,lwd=c(rep(4,4),NA,NA,6),cex=1.5)
			} else {
				legend('topleft',c('ATom 1 - Aug. 2016','Atom 2 - Feb. 2017','ATom 3 - Oct. 2017','ATom 4 - May 2018','Annual Mean'),col=c(cols,'black'),pch=c(rep(16,4),16),pt.cex=c(rep(1,4),1.5),pt.lwd=2,lwd=c(rep(4,4),6),cex=1.5)
			}
			dev.off()

			## only doing surface conc and gradients on annual means, so not plotting individual ATom surface or gradients here

		}
	}
	assign(paste(exp,'setcol',sep=''),expsetcol)
	assign(paste(exp,'meancol',sep=''),apply(expsetcol[,2:ncol(expsetcol)],1,mean))
	assign(paste(exp,'setgrad',sep=''),expsetgrad)
	assign(paste(exp,'meangrad',sep=''),apply(expsetgrad[,2:ncol(expsetgrad)],1,mean))
	assign(paste(exp,'setsurf',sep=''),expsetsurf)
	assign(paste(exp,'meansurf',sep=''),apply(expsetsurf[,2:ncol(expsetsurf)],1,mean))

	# make experiment plots with all models

	png(paste(modeldir,'/',exp,'_ATom_partial_columns.png',sep=''),width=1800,height=1200,pointsize=24)
	ylm=range(expsetcol[,2:ncol(expsetcol)],na.rm=T)
	plot(expsetcol$lat,expsetcol[,2],type='n',xlab='Latitude (deg N)',ylab='Delta CO2 (- MLO trend)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=paste(exp,' ATom 1000-300 mb 10-degree Partial Column Means',sep=''),ylim=ylm)
	for(mod in colnames(expsetcol)[2:ncol(expsetcol)]){
                points(expsetcol$lat,expsetcol[mod][,1],type='b',pch=16,cex=1,lwd=4,col=modcols[which(models==mod)])
	}
	points(expsetcol$lat,get(paste(exp,'meancol',sep='')),type='b',pch=16,cex=1.5,lwd=6,col='black')
        legend('topleft',c(colnames(expsetcol)[2:ncol(expsetcol)],'Mean'),col=c(modcols[which(is.element(models,colnames(expsetcol)[2:ncol(expsetcol)]))],'black'),pch=16,cex=1,lwd=c(rep(4,(ncol(expsetcol)-1)),6),pt.cex=c(rep(1,(ncol(expsetcol)-1)),1.5))
	dev.off()

	png(paste(modeldir,'/',exp,'_ATom_vertical_gradients.png',sep=''),width=1800,height=1200,pointsize=24)
        ylm=range(expsetgrad[,2:ncol(expsetgrad)],na.rm=T)
        plot(expsetgrad$lat,expsetgrad[,2],type='n',xlab='Latitude (deg N)',ylab='Delta CO2 (- MLO trend)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=paste(exp,' ATom 1000-',gradco,' minus ',gradco,'-300 mb 10-degree Vertical Gradients',sep=''),ylim=ylm)
        for(mod in colnames(expsetgrad)[2:ncol(expsetgrad)]){
                points(expsetgrad$lat,expsetgrad[mod][,1],type='b',pch=16,cex=1,lwd=4,col=modcols[which(models==mod)])
        }
        points(expsetgrad$lat,get(paste(exp,'meangrad',sep='')),type='b',pch=16,cex=1.5,lwd=6,col='black')
        legend('topleft',c(colnames(expsetgrad)[2:ncol(expsetgrad)],'Mean'),col=c(modcols[which(is.element(models,colnames(expsetcol)[2:ncol(expsetcol)]))],'black'),pch=16,cex=1,lwd=c(rep(4,(ncol(expsetgrad)-1)),6),pt.cex=c(rep(1,(ncol(expsetgrad)-1)),1.5))
        dev.off()

	png(paste(modeldir,'/',exp,'_ATom_surface_concentrations.png',sep=''),width=1800,height=1200,pointsize=24)
        ylm=range(expsetsurf[,2:ncol(expsetsurf)],na.rm=T)
        plot(expsetsurf$lat,expsetsurf[,2],type='n',xlab='Latitude (deg N)',ylab='Delta CO2 (- MLO trend)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=paste(exp,' ATom 1000-900 mb 10-degree Surface Concentrations',sep=''),ylim=ylm)
        for(mod in colnames(expsetsurf)[2:ncol(expsetsurf)]){
                points(expsetsurf$lat,expsetsurf[mod][,1],type='b',pch=16,cex=1,lwd=4,col=modcols[which(models==mod)])
        }
        points(expsetsurf$lat,get(paste(exp,'meansurf',sep='')),type='b',pch=16,cex=1.5,lwd=6,col='black')
        legend('topleft',c(colnames(expsetsurf)[2:ncol(expsetsurf)],'Mean'),col=c(modcols[which(is.element(models,colnames(expsetcol)[2:ncol(expsetcol)]))],'black'),pch=16,cex=1,lwd=c(rep(4,(ncol(expsetsurf)-1)),6),pt.cex=c(rep(1,(ncol(expsetsurf)-1)),1.5))
        dev.off()

}

# make annual mean plot with obs and experiment means

png(paste(modeldir,'/ATom_partial_columns_ann_w_exp_means.png',sep=''),width=1800,height=1200,pointsize=24)
ylm=range(c(ISmeancol,OGmeancol,LNLGmeancol,LNLGOGISmeancol,LNLGISmeancol),na.rm=T)
plot(ISsetcol$lat,ISmeancol,type='n',xlab='Latitude (deg N)',ylab='Delta CO2 (- MLO trend)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='ATom 1000-300 mb 10-degree Partial Column Means',ylim=ylm)
for(exp in experiments){
	points(expsetcol$lat,get(paste(exp,'meancol',sep='')),type='b',pch=16,cex=1,lwd=4,col=expcols[which(experiments==exp)])
}
points(obsanncol$lat,obsanncol$co2,type='b',pch=16,cex=1.5,lwd=6,col='black')
legend('topleft',c(experiments,'ATom Observation'),col=c(expcols,'black'),pch=16,lwd=c(rep(4,length(experiments)),6),pt.cex=c(rep(1,length(experiments)),1.5),cex=1.5)
dev.off()

# diff plot 
png(paste(modeldir,'/ATom_partial_columns_ann_exp_means_offsets.png',sep=''),width=1800,height=1200,pointsize=24)
ylm=c(-0.5,1.5)
plot(ISsetcol$lat,ISmeancol-obsanncol$co2,type='n',xlab='Latitude (deg N)',ylab='Delta CO2 (- ATom Obs)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='ATom 1000-300 mb 10-degree Partial Column Means',ylim=ylm)
abline(h=0)
for(exp in experiments){
	points(expsetcol$lat,get(paste(exp,'meancol',sep=''))-obsanncol$co2,type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(experiments==exp)])
}
legend('topleft',experiments,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()

# diff plot
png(paste(modeldir,'/ATom_partial_columns_ann_exp_means_offsets.png',sep=''),width=1800,height=1200,pointsize=24)
ylm=c(-0.5,1.5)
plot(ISsetcol$lat,ISmeancol-obsanncol$co2,type='n',xlab='Latitude (deg N)',ylab='Delta CO2 (- ATom Obs)',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main='ATom 1000-300 mb 10-degree Partial Column Means',ylim=ylm)
abline(h=0)
for(exp in experiments){
        points(expsetcol$lat,get(paste(exp,'meancol',sep=''))-obsanncol$co2,type='b',pch=16,cex=1.5,lwd=6,col=expcols[which(experiments==exp)])
}
legend('topleft',experiments,col=expcols,pch=16,lwd=6,pt.cex=1.5)
dev.off()

# redo for paper:

png(paste(modeldir,'/fig2b_colmeanoffsets.png',sep=''),width=650,height=500*12/15,pointsize=12) 
par(mar=c(4,5,4,1))
par(mgp=c(3,1,0))
#ylm=c(-0.5,1.5)
ylm=range(c(ISmeancol-obsanncol$co2,OGmeancol-obsanncol$co2,LNLGmeancol-obsanncol$co2,LNLGOGISmeancol-obsanncol$co2,LNLGISmeancol-obsanncol$co2),na.rm=T)
plot(ISsetcol$lat,ISmeancol-obsanncol$co2,type='n',xlab=expression(paste("Latitude (",degree, "N)")),ylab=expression(paste(Delta,CO[2],' (ppm)')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Column-mean model minus observation'),ylim=ylm)
abline(h=0)
for(exp in experiments){
        points(expsetcol$lat,get(paste(exp,'meancol',sep=''))-obsanncol$co2,type='b',pch=19,cex=1.5,lwd=4,col=expcols[which(experiments==exp)])
}
legend('topleft',experiments,col=expcols,pch=19,lwd=4,pt.cex=1.5,bty='n',inset=0.05)
dev.off()

png(paste(modeldir,'/fig2b_colmeanoffsets_windivmod.png',sep=''),width=650,height=500*12/15,pointsize=12) 
par(mar=c(4,5,4,1))
par(mgp=c(3,1,0))
#ylm=c(-0.5,1.5)
ylm=range(c(ISmeancol-obsanncol$co2,OGmeancol-obsanncol$co2,LNLGmeancol-obsanncol$co2,LNLGOGISmeancol-obsanncol$co2,LNLGISmeancol-obsanncol$co2),na.rm=T)
plot(ISsetcol$lat,ISmeancol-obsanncol$co2,type='n',xlab=expression(paste("Latitude (",degree, "N)")),ylab=expression(paste(Delta,CO[2],' (ppm)')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Column-mean model minus observation'),ylim=ylm)
abline(h=0)
for(exp in experiments){
	expsetcol=get(paste(exp,'setcol',sep=''))
	for(mod in colnames(expsetcol)[2:ncol(expsetcol)]){
                #lines(expsetcol$lat,expsetcol[mod][,1]-obsanncol$co2,lwd=1,col=modcols[which(models==mod)]) # diff cols by mod
                lines(expsetcol$lat,expsetcol[mod][,1]-obsanncol$co2,lwd=1,col=expcols[which(experiments==exp)]) # diff line types by mod?
	}
        points(expsetcol$lat,get(paste(exp,'meancol',sep=''))-obsanncol$co2,type='b',pch=19,cex=1.5,lwd=4,col=expcols[which(experiments==exp)])
}
legend('topleft',experiments,col=expcols,pch=19,lwd=4,pt.cex=1.5,bty='n',inset=0.05)
dev.off()

png(paste(modeldir,'/fig2b_colmeanoffsets_windivmod_shaded.png',sep=''),width=650,height=500*12/15,pointsize=12)
par(mar=c(4,5,4,1))
par(mgp=c(3,1,0))
#ylm=c(-0.5,1.5)
ylm=range(c(ISmeancol-obsanncol$co2,OGmeancol-obsanncol$co2,LNLGmeancol-obsanncol$co2,LNLGOGISmeancol-obsanncol$co2,LNLGISmeancol-obsanncol$co2),na.rm=T)
plot(ISsetcol$lat,ISmeancol-obsanncol$co2,type='n',xlab=expression(paste("Latitude (",degree, "N)")),ylab=expression(paste(Delta,CO[2],' (ppm)')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Column-mean model minus observation'),ylim=ylm)
abline(h=0)
for(exp in experiments){
        expsetcol=get(paste(exp,'setcol',sep=''))
	expsds=apply(expsetcol[,2:ncol(expsetcol)],1,sd)
        for(mod in colnames(expsetcol)[2:ncol(expsetcol)]){
		polygon(c(expsetcol$lat,rev(expsetcol$lat)),c(-1*expsds+get(paste(exp,'meancol',sep=''))-obsanncol$co2,rev(expsds+get(paste(exp,'meancol',sep=''))-obsanncol$co2)),border=NA,col=lighten(expcols[which(experiments==exp)],0.5))
        }
}
for(exp in experiments){
        points(expsetcol$lat,get(paste(exp,'meancol',sep=''))-obsanncol$co2,type='b',pch=19,cex=1.5,lwd=4,col=expcols[which(experiments==exp)])
}
legend('topleft',experiments,col=expcols,pch=19,lwd=4,pt.cex=1.5,bty='n',inset=0.05)
dev.off()

png(paste(modeldir,'/fig2b_colmeanoffsets_windivmod_sdrange.png',sep=''),width=650,height=500*12/15,pointsize=12)
par(mar=c(4,5,4,1))
par(mgp=c(3,1,0))
ylm=c(-0.35,1.3)
plot(ISsetcol$lat,ISmeancol-obsanncol$co2,type='n',xlab=expression(paste("Latitude (",degree, "N)")),ylab=expression(paste(CO[2],' Bias (ppm)')),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Column-mean model minus observation'),ylim=ylm)
abline(h=0)
for(exp in experiments){
        expsetcol=get(paste(exp,'setcol',sep=''))
        expsds=apply(expsetcol[,2:ncol(expsetcol)],1,sd)
        for(mod in colnames(expsetcol)[2:ncol(expsetcol)]){
                lines(expsetcol$lat,expsds+get(paste(exp,'meancol',sep=''))-obsanncol$co2,lwd=1,col=lighten(expcols[which(experiments==exp)],0.5),lty=3)
                lines(expsetcol$lat,-1*expsds+get(paste(exp,'meancol',sep=''))-obsanncol$co2,lwd=1,col=lighten(expcols[which(experiments==exp)],0.5),lty=3)
        }
}
for(exp in experiments){
        points(expsetcol$lat,get(paste(exp,'meancol',sep=''))-obsanncol$co2,type='b',pch=19,cex=1.5,lwd=4,col=expcols[which(experiments==exp)])
}
legend('topleft',experiments,col=expcols,pch=19,lwd=4,pt.cex=1.5,bty='n',inset=0.05)
dev.off()


## EPS

width_cm <- 8.9
width_in <- width_cm / 2.54
sc <- 0.21 # Scaling factor for absolute line widths
new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)
file_name <- paste(modeldir,'/fig2b_colmeanoffsets_windivmod_sdrange.eps',sep='')
cairo_ps(filename = file_name,
         width = width_in,
         height = width_in*500/650*12/15,
         pointsize = new_pointsize)

ylm=c(-0.35,1.3)
plot(ISsetcol$lat,ISmeancol-obsanncol$co2,type='n',xlab=expression(paste("Latitude (",degree, "N)")),ylab='',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,main=expression('b. Column-mean model minus observation'),ylim=ylm,axes=F)
box(lwd=0.5)
axis(1,cex.axis=1.5,lwd=0.5)
axis(2,cex.axis=1.5,lwd=0.5)
mtext(expression(paste(CO[2],' Bias (ppm)')),2,2.5,cex=1.5)
abline(h=0,lwd=0.5)
for(exp in experiments){
        expsetcol=get(paste(exp,'setcol',sep=''))
        expsds=apply(expsetcol[,2:ncol(expsetcol)],1,sd)
        for(mod in colnames(expsetcol)[2:ncol(expsetcol)]){
                lines(expsetcol$lat,expsds+get(paste(exp,'meancol',sep=''))-obsanncol$co2,lwd=0.5,col=expcols[which(experiments==exp)],lty=3)
                lines(expsetcol$lat,-1*expsds+get(paste(exp,'meancol',sep=''))-obsanncol$co2,lwd=0.5,col=expcols[which(experiments==exp)],lty=3)
        }
}
for(exp in experiments){
        points(expsetcol$lat,get(paste(exp,'meancol',sep=''))-obsanncol$co2,type='b',pch=19,cex=1.5,lwd=4*sc*1.5,col=expcols[which(experiments==exp)])
}
legend('topleft',experiments,col=expcols,pch=19,lwd=4*sc*1.5,pt.cex=1.5,bty='n',inset=0.05)
dev.off()

system(paste('cp ',modeldir,'/fig2b_colmeanoffsets_windivmod_sdrange.eps ',modeldir,'/fig2b.eps',sep=''))
