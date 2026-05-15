## Program to correlate model posterior CO2 and fluxes for each individual press x lat bin for 50 models against 50 mean fluxes for a selected region and plot results
##  also examines EOFs


library(abind)
library(pracma)
library('RColorBrewer')


## Options:

# define reference lat / prs for optional subtraction
reflat=19.4; refprs= 680 # corresponding to MLO 
refsub=F # whether to subtract global mean from each model
refsubopt='mlo' # subtract value interpolated to reflat and refprs from each model
#refsubopt='bg' # subtract global mean for each model

exclog=F # whether to exclude OG and LNLGOGIS
oneexponly=NULL # 'LNLG' # option to only do one experiment

rsqco=0 # option to filter on rsqd

latlim=c(-70,80) # trim concentrations
latlimfluxes=F # whether to trim fluxes to latlim

fitmodmeans=F; fitexpmeans=F # options instead of individual inversions, can not both be true

fillmissing=T # to interp any missing bins (none expected at 10x100 for in situ obs)

obsdir='OBS/BINFIT'

models=c('AMES','BAKER','CAMS','CMS-FLUX','COLA','CT','OU','TM5-4DVAR','UT','WOMBAT') 
experiments=c('IS','OG','LNLG','LNLGOGIS','LNLGIS') # ,'unopt') 
if(exclog) experiments=c('IS','LNLG','LNLGIS') # excluding OG
if(!is.null(oneexponly)) experiments=oneexponly
modeldir='V10MIP'

# set plot size and axis ranges
pngw=650;pngh=650;psize=12
xlim1=c(-70,80);ylim1=c(1000,300)
xtick=20
res='10x100';latspan=10


## Read in and stack 2D (P vs lat) annual mean fields for each model:

allmat=NULL # 2D fields expanded into single columns
allexparr=NULL # 2D fields of experiment means preserved (averaged for experiment)
modnames=NULL
expnames=NULL
allbg=NULL
latvals=scan(paste(modeldir,'/',models[1],'/',experiments[1],'/Annual_Mean_CO2_OP_',res,'.txt',sep=''),nlines=1)
prsvals=scan(paste(modeldir,'/',models[1],'/',experiments[1],'/Annual_Mean_CO2_OP_',res,'.txt',sep=''),skip=1,nlines=1)
latwt=cos(latvals/180*pi)
for(exp in experiments){
        print(exp)
	expmat=NULL # 2D fields expanded into single columns
	exparr=NULL # 2D fields preserved
	nmod=0

        for(mod in models){

		print(mod)
		nmod=nmod+1
		modnames=c(modnames,mod)
		expnames=c(expnames,exp)

		annmat=read.table(paste(modeldir,'/',mod,'/',exp,'/Annual_Mean_CO2_OP_',res,'.txt',sep=''),skip=2) # from parent directory binning by flt

		if(fillmissing&any(is.na(annmat))){ # interpolate by lat and then by prs to fill missing bins (e.g. for flasks)
			print(paste('Filling',sum(is.na(annmat)),'missing',mod,'/',exp,'bins'))
			for(i in c(1:ncol(annmat))){
				if(sum(!is.na(annmat[,i]))>1){
					annmat[,i][is.na(annmat[,i])]=approx(c(1:nrow(annmat))[!is.na(annmat[,i])],annmat[,i][!is.na(annmat[,i])],c(1:nrow(annmat))[is.na(annmat[,i])],rule=2)$y
				}
			}
			for(i in c(1:nrow(annmat))){
				if(sum(!is.na(annmat[i,]))>1){
					annmat[i,][is.na(annmat[i,])]=approx(c(1:ncol(annmat))[!is.na(annmat[i,])],annmat[i,][!is.na(annmat[i,])],c(1:ncol(annmat))[is.na(annmat[i,])],rule=2)$y
				}
			}
			print(paste(sum(is.na(annmat)),'remaining missing'))
		}

		allbg=c(allbg,sum(unlist(annmat*latwt))/sum(latwt*length(prsvals)))
		if(refsub){
			if(refsubopt=='mlo'){
				# subtract reference
				tmp=matrix(unlist(annmat),ncol=ncol(annmat),nrow=nrow(annmat)) # annmat is non-numeric list at this point
				refco2=interp2(prsvals,latvals,tmp,refprs,reflat,method="linear") # MLO
				annmat=annmat-refco2
			} else if(refsubopt=='bg'){
				# subtract mean
				annmat=annmat-sum(unlist(annmat*latwt))/sum(latwt*length(prsvals))
			}
		}

		if(is.null(allmat)){ # first model/exp combo
			allmat=unlist(annmat)
			alllat=rep(latvals,length(prsvals))
			allprs=rep(prsvals,each=length(latvals))
		} else {
			allmat=cbind(allmat,unlist(annmat))
		}
		if(is.null(expmat)){ # first model/exp combo for this experiment
			exparr=annmat
			expmat=unlist(annmat)
		} else {
			exparr=exparr+annmat
			expmat=expmat+unlist(annmat) 
		}

        }
	if(is.null(allexparr)){ # first exp
		allexparr=exparr/nmod
	} else {
		allexparr=abind(allexparr,exparr/nmod,along=3)
	}

}

# write out results for use by flux_cor_ecfc.r
write(c('latitude','pressure',paste(modnames,'_',expnames,sep='')),paste(modeldir,'/Annual_mean_2D_fields_as_columns.txt',sep=''),ncol=ncol(allmat)+2)
write(t(cbind(alllat,allprs,allmat)),paste(modeldir,'/Annual_mean_2D_fields_as_columns.txt',sep=''),ncol=ncol(allmat)+2,append=T)
write(allbg,paste(modeldir,'/Annual_mean_2D_field_averages.txt',sep=''),ncol=1)

# trim model concentrations for latlim
allmat=allmat[alllat>=latlim[1]&alllat<=latlim[2],]
allprs=allprs[alllat>=latlim[1]&alllat<=latlim[2]]
alllat=alllat[alllat>=latlim[1]&alllat<=latlim[2]]
allexparr=allexparr[latvals>=latlim[1]&latvals<=latlim[2],,]


## Do for obs too:

variables=c('CO2_OP','CO2_AO2_m_CO2_NOAA','CO2_QCLS_m_CO2_NOAA','CO2_MED_m_CO2_NOAA','CO2_PFP_m_CO2_NOAA')
obsannmat=NULL # 2D fields expanded into single columns
obsannarr=NULL # 2D fields preserved
obsannbg=NULL
for(var in variables){

        print(var)
        annmat=read.table(paste(obsdir,'/Annual_Mean_',var,'_',res,'.txt',sep=''),skip=2)

	if(fillmissing&any(is.na(annmat))){ # interpolate by lat and then by prs to fill missing bins (e.g. for flasks)
                print(paste('Filling',sum(is.na(annmat)),'missing',var,'bins'))
                for(i in c(1:ncol(annmat))){
                        if(sum(!is.na(annmat[,i]))>1){
                                annmat[,i][is.na(annmat[,i])]=approx(c(1:nrow(annmat))[!is.na(annmat[,i])],annmat[,i][!is.na(annmat[,i])],c(1:nrow(annmat))[is.na(annmat[,i])],rule=2)$y
                        }
                }
                for(i in c(1:nrow(annmat))){
                        if(sum(!is.na(annmat[i,]))>1){
                                annmat[i,][is.na(annmat[i,])]=approx(c(1:ncol(annmat))[!is.na(annmat[i,])],annmat[i,][!is.na(annmat[i,])],c(1:ncol(annmat))[is.na(annmat[i,])],rule=2)$y
                        }
                }
                print(paste(sum(is.na(annmat)),'remaining missing'))
        }

        obsannbg=c(obsannbg,sum(unlist(annmat*latwt))/sum(latwt*length(prsvals)))

        if(refsub){
                        if(refsubopt=='mlo'){
                                # subtract reference
                                tmp=matrix(unlist(annmat),ncol=ncol(annmat),nrow=nrow(annmat)) # annmat is non-numeric list at this point
                                refco2=interp2(prsvals,latvals,tmp,refprs,reflat,method="linear") # MLO
                                annmat=annmat-refco2
                        } else if(refsubopt=='bg'){
                                # subtract mean
                                annmat=annmat-sum(unlist(annmat*latwt))/sum(latwt*length(prsvals))
                        }
        }

        if(is.null(obsannmat)){ # first model/exp combo
                obsannarr=annmat
                obsannmat=unlist(annmat)
                obslat=rep(latvals,length(prsvals))
                obsprs=rep(prsvals,length(latvals))
        } else {
                obsannarr=abind(obsannarr,annmat,along=3)
                obsannmat=cbind(obsannmat,unlist(annmat))
        }

}

# write out results for use by flux_cor_ecfc.r
write(c('latitude','pressure',variables),paste(obsdir,'/Annual_mean_2D_fields_as_columns.txt',sep=''),ncol=ncol(obsannmat)+2)
write(t(cbind(obslat,obsprs,obsannmat)),paste(obsdir,'/Annual_mean_2D_fields_as_columns.txt',sep=''),ncol=ncol(obsannmat)+2,append=T)

# trim observed concentrations for latlim
obsannmat=obsannmat[obslat>=latlim[1]&obslat<=latlim[2],]
obslat=obslat[obslat>=latlim[1]&obslat<=latlim[2]]
# now trim latvals and prsvals too
latvals=latvals[latvals>=latlim[1]&latvals<=latlim[2]]


## Make model means:

modmeans=aggregate(t(allmat),by=list(mod=modnames),mean)
header=modmeans$mod
modmeans=t(modmeans[,2:ncol(modmeans)])
colnames(modmeans)=header
expcount=aggregate(t(!is.na(allmat)),by=list(mod=modnames),sum)
header=expcount$mod
expcount=t(expcount[,2:ncol(expcount)])
colnames(expcount)=header
# check if any are missing
if(any(expcount!=length(experiments))){
	print('Missing some model experiments in model means')
	stop()
}

# calculate experiment means (to exclude OU as allexparr includes OU OG in OG mean)
expmeans=aggregate(t(allmat),by=list(exp=expnames),mean)
header=expmeans$exp
if(!is.null(oneexponly)) expmeans=matrix(expmeans,nrow=length(experiments)) # to allow for running only one experiment
expmeans=expmeans[,2:ncol(expmeans)]
if(!is.null(oneexponly)) expmeans=matrix(expmeans,nrow=length(experiments)) # to allow for running only one experiment
expmeans=t(expmeans) # has to be a matrix to transpose
colnames(expmeans)=header
if(is.null(oneexponly)) expmeans=expmeans[,experiments] # put back in column order to match experiments


modcount=aggregate(t(!is.na(allmat)),by=list(exp=expnames),sum)
header=modcount$exp
modcount=t(modcount[,2:ncol(modcount)])
colnames(modcount)=header
modcount=modcount[,experiments] # put back in column order to match experiments
# check if any are missing
if(any(modcount!=length(models))){
        print('Missing some experiment models in experiment means')
        stop()
}


## Make plots:

# set plot colors
rdylblfunc=colorRampPalette(brewer.pal(11,'RdYlBu'))
rdylbl=rdylblfunc(127)
blylrd=rev(rdylbl)
usecol=blylrd # for calls to xsect_plot.r

rdwtblfunc=colorRampPalette(brewer.pal(11,'RdBu'))
rdwtbl=rdwtblfunc(127)
blwtrd=rev(rdwtbl)

cols=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2')) # 17 diff colors
cols[6]='gold' # yellow too hard to see
cols=c(cols,cols,cols) # now 51

expcols=cols[c(1,2,3,4,5)]
expcols=expcols[is.element(experiments,expnames)]

source('xsect_plot.r')

# plot experiment means
zlm=range(c(range(allexparr),range(obsannmat[,variables=='CO2_OP'])))
for(exp in experiments){ # loop and plot exp means

	png(paste(modeldir,'/Annual_Mean_',exp,'_average.png',sep=''),width=650,height=650,pointsize=12)
	print(exp)
	print(range(allexparr[,,which(experiments==exp)]))
	print(range(expmeans[,colnames(expmeans)==exp])) # identical ### column order of expmeans does not match experiments
	if(is.null(oneexponly)){
		xs(allexparr[,,which(experiments==exp)],color=blylrd,colbar=T,zlm=zlm,colbarint=F)
	} else {
		xs(as.matrix(allexparr),color=blylrd,colbar=T,zlm=zlm,colbarint=F)
	}
	dev.off()

}

# plot obs

png('Annual_Mean_Obs_average.png',width=650,height=650,pointsize=12) # this one has zlim matching experiment mean plots
print('Obs')
print(range(obsannmat[,variables=='CO2_OP'],na.rm=T))
xs(matrix(obsannmat[,variables=='CO2_OP'],ncol=length(prsvals),nrow=length(latvals)),color=blylrd,colbar=T,zlm=zlm,colbarint=F) 
dev.off()


# for paper (requires running plot_atom_col.r first)
zlm=range(obsannmat[,variables=='CO2_OP']) # this one zoomed in to just obs
png('Annual_Mean_Obs_average_fig2a_wcolbar.png',width=650,height=500,pointsize=12) # since colbar=T and colbarint=F gives heights=rep(c(12,2.0)
xs(matrix(obsannmat[,variables=='CO2_OP'],ncol=length(prsvals),nrow=length(latvals)),color=usecol,colbar=T,colbarint=F,zlm=zlm,main=expression('a. Observed'),legendarrows=F,units=expression(paste(Delta,CO[2],' (ppm)')))
dev.off()
system('montage Annual_Mean_Obs_average_fig2a_wcolbar.png V10MIP/fig2b_colmeanoffsets_windivmod_sdrange.png -geometry +1+3 -gravity North fig2_annxsect_colmeanoffsets.png')

## EPS
width_cm <- 8.9
width_in <- width_cm / 2.54
sc <- 0.21 # Scaling factor for absolute line widths
new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)
file_name <- 'Annual_Mean_Obs_average_fig2a_wcolbar.eps'
cairo_ps(filename = file_name,
         width = width_in,
         height = width_in*500/650,
         pointsize = new_pointsize)
xs(matrix(obsannmat[,variables=='CO2_OP'],ncol=length(prsvals),nrow=length(latvals)),color=usecol,colbar=T,colbarint=F,zlm=zlm,main=expression('a. Observed'),legendarrows=F,units=expression(paste(Delta,CO[2],' (ppm)')),axlwd=0.5)
dev.off()
system('cp Annual_Mean_Obs_average_fig2a_wcolbar.eps fig2a.eps')


png(paste(modeldir,'/All_model_annmns_1-5.png',sep=''),width=pngw*4,height=pngh*3,pointsize=psize)
par(mfrow=c(5,5))
colorder=c(which(expnames=='OG'),which(expnames=='LNLGOGIS'),which(expnames=='IS'),which(expnames=='LNLGIS'),which(expnames=='LNLG'))
for(j in colorder[c(1:5,11:15,21:25,31:35,41:45)]){
        plotvar=allmat[,j]
        plotvar[plotvar>4]=4
        plotvar[plotvar<(-4)]=-4
        xs(matrix(plotvar,ncol=length(prsvals),nrow=length(latvals)),zlm=c(-4,4),main=paste(expnames[j],modnames[j]),color=blylrd,colbar=F,txtcex=2)
}
dev.off()

png(paste(modeldir,'/All_model_annmns_6-10.png',sep=''),width=pngw*4,height=pngh*3,pointsize=psize)
par(mfrow=c(5,5))
colorder=c(which(expnames=='OG'),which(expnames=='LNLGOGIS'),which(expnames=='IS'),which(expnames=='LNLGIS'),which(expnames=='LNLG'))
for(j in colorder[c(6:10,16:20,26:30,36:40,46:50)]){
        plotvar=allmat[,j]
        plotvar[plotvar>4]=4
        plotvar[plotvar<(-4)]=-4
        xs(matrix(plotvar,ncol=length(prsvals),nrow=length(latvals)),zlm=c(-4,4),main=paste(expnames[j],modnames[j]),color=blylrd,colbar=F,txtcex=2)
}
dev.off()

png(paste(modeldir,'/All_model_minus_obs_annmns_1-5.png',sep=''),width=pngw*4,height=pngh*3,pointsize=psize)
par(mfrow=c(5,5))
colorder=c(which(expnames=='OG'),which(expnames=='LNLGOGIS'),which(expnames=='IS'),which(expnames=='LNLGIS'),which(expnames=='LNLG'))
for(j in colorder[c(1:5,11:15,21:25,31:35,41:45)]){
        plotvar=allmat[,j]-obsannmat[,variables=='CO2_OP']
        plotvar[plotvar>2]=2
        plotvar[plotvar<(-2)]=-2
        xs(matrix(plotvar,ncol=length(prsvals),nrow=length(latvals)),zlm=c(-2,2),main=paste(expnames[j],modnames[j]),color=blwtrd,colbar=F)
}
dev.off()

png(paste(modeldir,'/All_model_minus_obs_annmns_6-10.png',sep=''),width=pngw*4,height=pngh*3,pointsize=psize)
par(mfrow=c(5,5))
colorder=c(which(expnames=='OG'),which(expnames=='LNLGOGIS'),which(expnames=='IS'),which(expnames=='LNLGIS'),which(expnames=='LNLG'))
for(j in colorder[c(6:10,16:20,26:30,36:40,46:50)]){
        plotvar=allmat[,j]-obsannmat[,variables=='CO2_OP']
        plotvar[plotvar>2]=2
        plotvar[plotvar<(-2)]=-2
        xs(matrix(plotvar,ncol=length(prsvals),nrow=length(latvals)),zlm=c(-2,2),main=paste(expnames[j],modnames[j]),color=blwtrd,colbar=F)
}
dev.off()


png(paste(modeldir,'/All_model_minus_obs_annmns_colorbar.png',sep=''),height=200,width=500,pointsize=24)

zlm=c(-2,2)
xl1=zlm[1];xl2=zlm[2]; XX=seq(xl1,xl2,length.out=length(blwtrd)); # XX is edges of cells, XMAT is center
XMAT=t(matrix(XX[2:length(XX)]-(XX[2]-XX[1])/2,nrow=1,byrow=T) ) #XMAT is the matrix used for the image of the color bar
par(mar=c(3,1,1,1)) # impact of margins depends on height and width and psize of main plot
par(mgp=c(1.7,1,0))
squashhz=0.2
squashvt=1.2
image(x=XX,y=c(0,1),z=XMAT,zlim=zlm,xlim=c(zlm[1]-squashhz*diff(zlm),zlm[2]+squashhz*diff(zlm)),ylim=c(-squashvt,1+squashvt),axes=F,ylab="",col=blwtrd,xlab=expression(paste(Delta,CO[2],' (ppm)')),cex.lab=1.5)

legendarrows=T
if(legendarrows){
        arrowwid=(xl2-xl1)*0.05
        polygon(c(xl2,xl2+arrowwid,xl2,xl2),c(0,0.5,1,0),col=tail(blwtrd,1),border=NA)
        polygon(c(xl1,xl1-arrowwid,xl1,xl1),c(0,0.5,1,0),col=head(blwtrd,1),border=NA)
        lines(c(xl2,xl2+arrowwid,xl2),c(0,0.5,1))
        lines(c(xl1,xl1-arrowwid,xl1),c(0,0.5,1))
        lines(c(xl1,xl2),c(0,0))
        lines(c(xl1,xl2),c(1,1))
} else {
        lines(c(xl1,xl2,xl2,xl1,xl1),c(1,1,0,0,1))
}

breaks.axis=NULL # can specify color bar tick locations
if(is.null(breaks.axis)){ aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)] }else{AA=breaks.axis}
axis(side=1,pos=0,at=AA,cex.axis=1.5)

dev.off()

png(paste(modeldir,'/All_model_annmns_colorbar.png',sep=''),height=200,width=500,pointsize=24)

zlm=c(-4,4)
xl1=zlm[1];xl2=zlm[2]; XX=seq(xl1,xl2,length.out=length(blylrd)); # XX is edges of cells, XMAT is center
XMAT=t(matrix(XX[2:length(XX)]-(XX[2]-XX[1])/2,nrow=1,byrow=T) ) #XMAT is the matrix used for the image of the color bar
par(mar=c(3,1,1,1)) # impact of margins depends on height and width and psize of main plot
par(mgp=c(1.7,1,0))
squashhz=0.2
squashvt=1.2
image(x=XX,y=c(0,1),z=XMAT,zlim=zlm,xlim=c(zlm[1]-squashhz*diff(zlm),zlm[2]+squashhz*diff(zlm)),ylim=c(-squashvt,1+squashvt),axes=F,ylab="",col=blylrd,xlab=expression(paste(Delta,CO[2],' (ppm)')),cex.lab=1.5)

legendarrows=T
if(legendarrows){
        arrowwid=(xl2-xl1)*0.05
        polygon(c(xl2,xl2+arrowwid,xl2,xl2),c(0,0.5,1,0),col=tail(blylrd,1),border=NA)
        polygon(c(xl1,xl1-arrowwid,xl1,xl1),c(0,0.5,1,0),col=head(blylrd,1),border=NA)
        lines(c(xl2,xl2+arrowwid,xl2),c(0,0.5,1))
        lines(c(xl1,xl1-arrowwid,xl1),c(0,0.5,1))
        lines(c(xl1,xl2),c(0,0))
        lines(c(xl1,xl2),c(1,1))
} else {
        lines(c(xl1,xl2,xl2,xl1,xl1),c(1,1,0,0,1))
}

breaks.axis=NULL # can specify color bar tick locations
if(is.null(breaks.axis)){ aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)] }else{AA=breaks.axis}
axis(side=1,pos=0,at=AA,cex.axis=1.5)

dev.off()


if(any(experiments=='OG')&any(experiments=='LNLG')){
	ogminlnlg=expmeans[,colnames(expmeans)=='OG']-expmeans[,colnames(expmeans)=='LNLG']

	png(paste(modeldir,'/Annual_Mean_OG_minus_LNLG_average.png',sep=''),width=650,height=650,pointsize=12)
	xs(matrix(ogminlnlg,ncol=length(prsvals),nrow=length(latvals)),zlm=range(ogminlnlg),color=blylrd,colbar=T,colbarint=F)
	dev.off()

}

# make model mean of obs-mod diffs
modmeanobsdiffs=aggregate(t(allmat-obsannmat[,variables=='CO2_OP']),by=list(mod=modnames),mean)
header=modmeanobsdiffs$mod
modmeanobsdiffs=t(modmeanobsdiffs[,2:ncol(modmeanobsdiffs)])
colnames(modmeanobsdiffs)=header

png(paste(modeldir,'/All_model_mean_conc_obsdiffs.png',sep=''),width=pngw*4,height=pngh*3,pointsize=psize)
par(mfrow=c(3,4))
for(j in c(1:ncol(modmeanobsdiffs))){
        xs(matrix(modmeanobsdiffs[,j],ncol=length(prsvals),nrow=length(latvals)),zlm=c(-1.3,1.25),main=paste(colnames(modmeanobsdiffs)[j]))
}
dev.off()


## Calculate EOFs:

scalefactor=rep(cos(latvals/180*pi),length(prsvals)) # converting ppm to PgC equiv with cos(lat)

pca=prcomp(t(allmat),scale.=1/scalefactor) # Proportion of Variance 0.4746 0.1517 0.08592 0.05853

eof1=matrix(pca$rotation[,1],ncol=length(prsvals),nrow=length(latvals))
eof2=matrix(pca$rotation[,2],ncol=length(prsvals),nrow=length(latvals))
eof3=matrix(pca$rotation[,3],ncol=length(prsvals),nrow=length(latvals))
eof4=matrix(pca$rotation[,4],ncol=length(prsvals),nrow=length(latvals))

# plot Proportion of Variance
png(paste(modeldir,'/EOF_Proportion_of_Variance_for_Conc.png',sep=''),width=pngw,height=pngh,pointsize=psize)
plot(c(1:min(dim(allmat))),summary(pca)$importance[2,],main='Proportion of Variance',xlab='EOF#',ylab='%')
print(paste('Conc EOF1 explains',summary(pca)$importance[2,1]))
dev.off()

# plot first 4 EOFs
png(paste(modeldir,'/EOFs_of_Conc_1-4.png',sep=''),width=pngw*2,height=pngh*2,pointsize=psize)
par(mfrow=c(2,2))
for(i in c(1:4)){
	xs(get(paste('eof',i,sep='')),zlm=c(-1,1)*max(abs(get(paste('eof',i,sep='')))),maintxt=paste('Empirical Orthogonal Function ',i,sep=''),color=blylrd)
}
dev.off()


## Calculate gridded SDs across models:

sdmat=matrix(apply(allmat,1,sd),ncol=length(prsvals),nrow=length(latvals))

# plot gridded SDs
png(paste(modeldir,'/Cross_Model_SDs.png',sep=''),width=pngw,height=pngh,pointsize=psize)
xs(sdmat,maintxt=paste('Standard Deviations across Inversions',sep=''),zlm=range(c(0,as.vector(sdmat)),na.rm=T),colbarint=F)
dev.off()

# plot gridded EOF1
png(paste(modeldir,'/Cross_Model_EOF1.png',sep=''),width=pngw,height=pngh,pointsize=psize)
xs(eof1,maintxt=paste('EOF1 across Inversions',sep=''),zlm=range(c(0,as.vector(eof1)),na.rm=T),colbarint=F,units='')
dev.off()


## Read in fluxes output by V10MIP/flux_v_lat.r:

fluxes=read.table(paste(modeldir,'/oco2_v10mip_total_fluxes_bylat_201606-201805.txt',sep=''),header=T,stringsAsFactors=F)

## for MIROC4-ACTM, need to match conc to fluxes (and exclude any conc without fluxes)
sel=is.element(paste(modnames,'_',expnames,sep=''),paste(fluxes[,2],'_',fluxes[,1],sep='')) # trim MIROC first
allmat=allmat[,sel]
modnames=modnames[sel]
expnames=expnames[sel]
allbg=allbg[sel]
fluxes=fluxes[is.element(paste(fluxes[,2],'_',fluxes[,1],sep=''),paste(modnames,'_',expnames,sep='')),]
# fix/ensure matching order
mch=match(paste(modnames,'_',expnames,sep=''),paste(fluxes[,2],'_',fluxes[,1],sep='')) # trim fluxes if needed second
allmat=allmat[,mch]
modnames=modnames[mch]
expnames=expnames[mch]

if(any(modnames!=fluxes[,2])){ stop('MISMATCH BETWEEN CONC AND FLUX ROWS') }
if(any(expnames!=fluxes[,1])){ stop('MISMATCH BETWEEN CONC AND FLUX ROWS') }

fluxlats=substr(names(fluxes),nchar(names(fluxes))-2,nchar(names(fluxes)))
fluxlats=fluxlats[3:length(fluxlats)]
fluxlats=gsub('\\.','-',fluxlats)
fluxlats=gsub('X','',fluxlats)
fluxlats=as.numeric(fluxlats)

fluxnums=t(fluxes[,3:ncol(fluxes)]) 

## do flux PCA 
# flux already in PgC/yr so does not need lat scaling
fluxpca=prcomp(t(fluxnums)) # returns 18 in summary(fluxpca)$importance[2,] corresponding to min(n,p)


## plot first four EOFs
feof1=fluxpca$rotation[,1]
feof2=fluxpca$rotation[,2]
feof3=fluxpca$rotation[,3]
feof4=fluxpca$rotation[,4]

ylm=range(c(fluxpca$rotation[,1:4]))
png(paste(modeldir,'/EOFs_of_Total_Fluxes.png',sep=''),width=pngw,height=pngh,pointsize=psize)
plot(fluxlats,feof1,type='b',main='Flux EOFs 1-4',xlab='Latitude',ylab='PgC/yr',col='red',ylim=ylm,lwd=2)
lines(fluxlats,feof2,col='blue',type='b',lwd=2)
lines(fluxlats,feof3,col='purple',type='b',lwd=2)
lines(fluxlats,feof4,col='orange',type='b',lwd=2)
legend('topleft',c('EOF1','EOF2','EOF3','EOF4'),col=c('red','blue','purple','orange'),pch=1,lwd=2)
dev.off()

png(paste(modeldir,'/EOF1_of_Total_Fluxes.png',sep=''),width=pngw,height=pngh/2,pointsize=psize*1.5)
par(mgp=c(2.5,1,0))
plot(fluxlats[fluxlats>=latlim[1]&fluxlats<=latlim[2]],feof1[fluxlats>=latlim[1]&fluxlats<=latlim[2]],type='b',main='OCO-2 v10 MIP Total Flux EOF1',xlab='Latitude',ylab='',col='red',lwd=3,pch=16,bg='red',cex=1.5)
abline(h=0)
dev.off()

# plot flux SDs
fluxsds=apply(fluxnums,1,sd)
png(paste(modeldir,'/Cross_Model_Flux_SDs.png',sep=''),width=pngw,height=pngh/2,pointsize=psize*1.5)
par(mgp=c(2,0.75,0))
plot(fluxlats,fluxsds,type='b',main='OCO-2 v10 MIP Total Flux SDs',xlab='Latitude',ylab='PgC/yr',col='red',pch=16,bg='red',cex=1.5,lwd=3)
dev.off()

# plot SD and EOF1
png(paste(modeldir,'/figS3_sdandeof1_cd.png',sep=''),width=650*2,height=500/2,pointsize=12)
par(mar=c(4,5,4,1))
par(mfrow=c(1,2))
plot(fluxlats[fluxlats>=latlim[1]&fluxlats<=latlim[2]],fluxsds[fluxlats>=latlim[1]&fluxlats<=latlim[2]],type='b',main=expression('c. Total Flux SDs'),xlab=expression(paste("Latitude (",degree, "N)")),ylab=expression(paste('PgC ',yr^-1)),col=cols[14],lwd=3,pch=16,bg=cols[14],cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
plot(fluxlats[fluxlats>=latlim[1]&fluxlats<=latlim[2]],feof1[fluxlats>=latlim[1]&fluxlats<=latlim[2]],main=expression('d. Total Flux EOF1'),xlab=expression(paste("Latitude (",degree, "N)")),ylab='',type='n',cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
abline(h=0) 
points(fluxlats[fluxlats>=latlim[1]&fluxlats<=latlim[2]],feof1[fluxlats>=latlim[1]&fluxlats<=latlim[2]],col=cols[14],lwd=3,pch=16,bg=cols[14],cex=1.5,type='b') # 3, 10, 14 are greens
dev.off()

# plot Proportion of Variance
png(paste(modeldir,'/EOF_Proportion_of_Variance_for_Fluxes.png',sep=''),width=pngw,height=pngh,pointsize=psize)
plot(c(1:min(dim(fluxnums))),summary(fluxpca)$importance[2,],main='Proportion of Variance',xlab='EOF#',ylab='%')
print(paste('Flux EOF1 explains',summary(fluxpca)$importance[2,1]))
dev.off()

ylm=range(fluxnums)
png(paste(modeldir,'/Cross_Model_Fluxes.png',sep=''),width=pngw,height=pngh/2,pointsize=psize*1.5)
par(mgp=c(2,0.75,0))
plot(fluxlats,fluxnums[,1],type='n',main='OCO-2 v10 MIP Total Fluxes',xlab='Latitude',ylab='PgC/yr',col='red',ylim=ylm,lwd=1.5)
for(i in c(1:ncol(fluxnums))){
	points(fluxlats,fluxnums[,i],type='l',col=cols[i],pch=16,bg=i,cex=1,lwd=3)
}
dev.off()

if(latlimfluxes){

	# binfit was done from 70S to 80N but fluxes go from 90S to 90N
	fluxnumsinclnp=fluxnums[fluxlats>=min(latvals),]
	fluxnums=fluxnums[fluxlats>=min(latvals)&fluxlats<=max(latvals),]
	fluxlats=fluxlats[fluxlats>=min(latvals)&fluxlats<=max(latvals)]

}


# 6 regions = SET, T, NET, SET-T, T-(SET+NET), NET-T
setflux=apply(fluxnums[fluxlats<=(-20),],2,sum)
tropflux=apply(fluxnums[fluxlats>=(-20)&fluxlats<=20,],2,sum)
netflux=apply(fluxnums[fluxlats>=20,],2,sum)
setmtropflux=setflux-tropflux
tropmsetnetflux=tropflux-(setflux+netflux)
netmtropflux=netflux-tropflux
regfluxes=cbind(setflux,tropflux,netflux,setmtropflux,tropmsetnetflux,netmtropflux)
labels=c('SET','Tropical','NET','SET-Trop','Trop-(SET+NET)','NET-Trop')
labels2=format(labels,width=max(nchar(labels)))

# make model mean regfluxes
modmeanregfluxes=aggregate(regfluxes,by=list(mod=modnames),mean)
modmeanmods=modmeanregfluxes[,1]
modmeanregfluxes=modmeanregfluxes[,2:ncol(modmeanregfluxes)]
expcountregfluxes=aggregate(!is.na(regfluxes),by=list(mod=modnames),sum)
expcountregfluxes=expcountregfluxes[,2:ncol(expcountregfluxes)]
# check if any others are missing
if(any(expcountregfluxes!=length(experiments))){
        print('Missing some model experiments in model mean reg fluxes')
        stop()
}

# check correlation to RMS
## better vs glb rather than vs individual region (with the exception of net which is slightly worse)
## better to unweighted rmsdev (with the exception of net which is worse)
rmsdevglb=NULL # ; rmsdevnet=NULL; rmsdevtrp=NULL; rmsdevset=NULL
rmsdevnet=NULL # ; rmsdevnet=NULL; rmsdevtrp=NULL; rmsdevset=NULL
rmsdevtrp=NULL # ; rmsdevnet=NULL; rmsdevtrp=NULL; rmsdevset=NULL
rmsdevset=NULL # ; rmsdevnet=NULL; rmsdevtrp=NULL; rmsdevset=NULL
for(j in c(1:length(expnames))){
        modmobs=allmat[,j]-obsannmat[,variables=='CO2_OP']
	rmsdevglb=c(rmsdevglb,sqrt(mean(modmobs^2)))
	rmsdevnet=c(rmsdevnet,sqrt(mean(modmobs[alllat>20]^2)))
	rmsdevtrp=c(rmsdevtrp,sqrt(mean(modmobs[alllat<20&alllat>(-20)]^2)))
	rmsdevset=c(rmsdevset,sqrt(mean(modmobs[alllat<(-20)]^2)))
}

png(paste(modeldir,'/NET_vs_rmsdevnet.png',sep=''),height=pngh,width=pngw,pointsize=psize)
par(mfrow=c(1,1)) # length(plotsel)))
par(mgp=c(2.75,1,0))
par(mar=c(5,4.5,2,1.0)+0.1)
par(oma=c(0,0,3,0))

plot(rmsdevnet,netflux,main='NET Flux vs. NET RMS Deviation',type='n',ylab=expression(paste('PgC ',yr^-1)),xlab='RMS Dev (ppm)',cex.main=2,cex.axis=1.5,cex.lab=1.5)
for(mod in models){
	points(rmsdevnet[modnames==mod],netflux[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
}
legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
legend('bottomright',experiments,col=expcols[c(1:length(experiments))],pt.cex=2.0,pt.lwd=3,pch=rep(16,length(experiments)),pt.bg='white',bty='n',lwd=rep(F,length(experiments)),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
dev.off()

png(paste(modeldir,'/SET_vs_rmsdevset.png',sep=''),height=pngh,width=pngw,pointsize=psize)
par(mfrow=c(1,1)) # length(plotsel)))
par(mgp=c(2.75,1,0))
par(mar=c(5,4.5,2,1.0)+0.1)
par(oma=c(0,0,3,0))
plot(rmsdevset,setflux,main='SET Flux vs. SET RMS Deviation',type='n',ylab=expression(paste('PgC ',yr^-1)),xlab='RMS Dev (ppm)',cex.main=2,cex.axis=1.5,cex.lab=1.5)
for(mod in models){
        points(rmsdevset[modnames==mod],setflux[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
}
legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
legend('bottomright',experiments,col=expcols[c(1:length(experiments))],pt.cex=2.0,pt.lwd=3,pch=rep(16,length(experiments)),pt.bg='white',bty='n',lwd=rep(F,length(experiments)),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
dev.off()

png(paste(modeldir,'/TROP_vs_rmsdevtrp.png',sep=''),height=pngh,width=pngw,pointsize=psize)
par(mfrow=c(1,1)) # length(plotsel)))
par(mgp=c(2.75,1,0))
par(mar=c(5,4.5,2,1.0)+0.1)
par(oma=c(0,0,3,0))
plot(rmsdevtrp,tropflux,main='TROP Flux vs. TROP RMS Deviation',type='n',ylab=expression(paste('PgC ',yr^-1)),xlab='RMS Dev (ppm)',cex.main=2,cex.axis=1.5,cex.lab=1.5)
for(mod in models){
        points(rmsdevtrp[modnames==mod],tropflux[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
}
legend('bottomleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
legend('topright',experiments,col=expcols[c(1:length(experiments))],pt.cex=2.0,pt.lwd=3,pch=rep(16,length(experiments)),pt.bg='white',bty='n',lwd=rep(F,length(experiments)),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
dev.off()

png(paste(modeldir,'/NET_vs_rmsdevglb.png',sep=''),height=pngh,width=pngw,pointsize=psize)
par(mfrow=c(1,1)) # length(plotsel)))
par(mgp=c(2.75,1,0))
par(mar=c(5,4.5,2,1.0)+0.1)
par(oma=c(0,0,3,0))

plot(rmsdevglb,netflux,main='NET Flux vs. Global RMS Deviation',type='n',ylab=expression(paste('PgC ',yr^-1)),xlab='RMS Dev (ppm)',cex.main=2,cex.axis=1.5,cex.lab=1.5)
for(mod in models){
        points(rmsdevglb[modnames==mod],netflux[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
}
legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
legend('bottomright',experiments,col=expcols[c(1:length(experiments))],pt.cex=2.0,pt.lwd=3,pch=rep(16,length(experiments)),pt.bg='white',bty='n',lwd=rep(F,length(experiments)),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
dev.off()

png(paste(modeldir,'/SET_vs_rmsdevglb.png',sep=''),height=pngh,width=pngw,pointsize=psize)
par(mfrow=c(1,1)) # length(plotsel)))
par(mgp=c(2.75,1,0))
par(mar=c(5,4.5,2,1.0)+0.1)
par(oma=c(0,0,3,0))
plot(rmsdevglb,setflux,main='SET Flux vs. Global RMS Deviation',type='n',ylab=expression(paste('PgC ',yr^-1)),xlab='RMS Dev (ppm)',cex.main=2,cex.axis=1.5,cex.lab=1.5)
for(mod in models){
        points(rmsdevglb[modnames==mod],setflux[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
}
legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
legend('bottomright',experiments,col=expcols[c(1:length(experiments))],pt.cex=2.0,pt.lwd=3,pch=rep(16,length(experiments)),pt.bg='white',bty='n',lwd=rep(F,length(experiments)),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
dev.off()

png(paste(modeldir,'/TROP_vs_rmsdevglb.png',sep=''),height=pngh,width=pngw,pointsize=psize)
par(mfrow=c(1,1)) # length(plotsel)))
par(mgp=c(2.75,1,0))
par(mar=c(5,4.5,2,1.0)+0.1)
par(oma=c(0,0,3,0))
plot(rmsdevglb,tropflux,main='TROP Flux vs. Global RMS Deviation',type='n',ylab=expression(paste('PgC ',yr^-1)),xlab='RMS Dev (ppm)',cex.main=2,cex.axis=1.5,cex.lab=1.5)
for(mod in models){
        points(rmsdevglb[modnames==mod],tropflux[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
}
legend('bottomleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
legend('topright',experiments,col=expcols[c(1:length(experiments))],pt.cex=2.0,pt.lwd=3,pch=rep(16,length(experiments)),pt.bg='white',bty='n',lwd=rep(F,length(experiments)),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
dev.off()




# make experiment mean regfluxes
tmp=regfluxes
modnamesel=rep(T,length(modnames))
tmp=tmp[modnamesel,]
expmeanregfluxes=aggregate(tmp,by=list(exp=expnames[modnamesel]),mean)
expmeanexps=expmeanregfluxes[,1]
expmeanregfluxes=expmeanregfluxes[,2:ncol(expmeanregfluxes)]
expmeanregfluxes=expmeanregfluxes[match(experiments,expmeanexps),] # put back in order to match experiments
modcountregfluxes=aggregate(!is.na(tmp),by=list(exp=expnames[modnamesel]),sum)
modcountexps=modcountregfluxes[,1]
modcountregfluxes=modcountregfluxes[,2:ncol(expcountregfluxes)]
modcountregfluxes=modcountregfluxes[match(experiments,modcountexps),] # put back in order to match experiments

if(res=='10x100'){ # would have to change prssel and latsel below to get to work for 5x50

	## plot example fit for a single bin
	pngdim=1500
	pngpts=30
	png(paste(modeldir,'/Cross_Model_Example_Fit.png',sep=''),height=pngdim,width=pngdim,pointsize=pngpts)
	par(mgp=c(2.5,1,0))
	par(mgp=c(2,0.75,0))
	par(mar=c(5,4,2,2)+0.1)
	par(oma=c(0,0,3,0))

	prssel=750
	latsel=35
	xvar=netflux; yvar=allmat[alllat==latsel&allprs==prssel,]
	xlm=range(xvar); ylm=range(yvar)

	plot(xvar,yvar,main='By Experiment',type='n',ylab='Annual Mean Flux (PgC/yr)',xlab='Concentration Metric (ppm)',cex.main=2,cex.axis=1.5,cex.lab=1.5,xlim=xlm,ylim=ylm)
	y1=par('usr')[3]; y2=par('usr')[4]
	for(mod in models){
		points(xvar[modnames==mod],yvar[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
	}
	abline(lm(yvar ~ xvar),col='Black',lwd=6)
	expmeanconc=aggregate(xvar,by=list(exp=expnames),mean)$x
	expsdconc=aggregate(xvar,by=list(exp=expnames),sd)$x
	expmeanflux=aggregate(yvar,by=list(exp=expnames),mean)$x
	expsdflux=aggregate(yvar,by=list(exp=expnames),sd)$x
	points(expmeanconc,expmeanflux,cex=2,col=expcols[1:length(experiments)],pch=16)
	segments(expmeanconc-expsdconc,expmeanflux,expmeanconc+expsdconc,expmeanflux,lwd=4,col=expcols[1:length(experiments)])
	segments(expmeanconc,expmeanflux-expsdflux,expmeanconc,expmeanflux+expsdflux,lwd=4,col=expcols[1:length(experiments)])
	legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=2,bg='white',box.col='white')
	legend('bottomright',experiments,col=expcols[c(1:length(experiments))],pt.cex=2.0,pt.lwd=3,pch=rep(16,length(experiments)),pt.bg='white',bty='n',lwd=rep(F,length(experiments)),bg='white',box.col='white')
	mtext(paste('SD =',round(sd(yvar),2)),3,-2)
	mtext(paste('res SD =',round(sd(lm(yvar ~ xvar)$resid),2)),3,-4)

	dev.off()

} # if(res=='10x100')


## Loop on lat+prs and region and do correlations:

rsqds=matrix(NA,ncol=6,nrow=nrow(allmat))
rs=matrix(NA,ncol=6,nrow=nrow(allmat))
icpts=matrix(NA,ncol=6,nrow=nrow(allmat))
coefs=matrix(NA,ncol=6,nrow=nrow(allmat))
invcoefs=matrix(NA,ncol=6,nrow=nrow(allmat))
ressd=matrix(NA,ncol=6,nrow=nrow(allmat))
invressd=matrix(NA,ncol=6,nrow=nrow(allmat))
uncred=matrix(NA,ncol=6,nrow=nrow(allmat))
for(i in c(1:nrow(allmat))){ # loop on conc bins
	for(j in c(1:6)){
		if(fitmodmeans){
			fit=lm(modmeans[i,] ~ modmeanregfluxes[,j])
			invfit=lm(modmeanregfluxes[,j] ~ modmeans[i,])
			rs[i,j]=cor(modmeans[i,],modmeanregfluxes[,j])
		} else if(fitexpmeans){
			fit=lm(expmeans[i,] ~ expmeanregfluxes[,j])
			invfit=lm(expmeanregfluxes[,j] ~ expmeans[i,])
			rs[i,j]=cor(expmeans[i,],expmeanregfluxes[,j])
		} else {
			fit=lm(allmat[i,] ~ regfluxes[,j])
			invfit=lm(regfluxes[,j] ~ allmat[i,])
			rs[i,j]=cor(allmat[i,],regfluxes[,j])
		}
		rsqds[i,j]=summary(fit)$r.squared
		icpts[i,j]=fit$coef[1]
		coefs[i,j]=fit$coef[2]
		invcoefs[i,j]=invfit$coef[2]
		ressd[i,j]=sd(fit$resid)
		invressd[i,j]=sd(invfit$resid)
		uncred[i,j]=(1-sd(fit$resid)/sd(allmat[i,]))*100
	}
}

# plot rsqds
png(paste(modeldir,'/Cross_Model_R-squared_Conc_vs_Flux_byreg.png',sep=''),width=pngw*3,height=pngh*2,pointsize=psize)

tmprsqds=rsqds
excllor2=F
if(excllor2){ # only plot cells with R^2>=rsqco
	tmprsqds[rsqds<rsqco]=NA
} 

par(mfrow=c(2,3))
for(i in c(1:6)){
	xs(matrix(tmprsqds[,i],ncol=length(prsvals),nrow=length(latvals)),zlm=range(rsqds),maintxt=labels[i],units='r-squared',lor2dots=T)
}
dev.off()


png(paste(modeldir,'/Cross_Model_R_Conc_vs_Flux_byreg.png',sep=''),width=pngw*3,height=pngh*2,pointsize=psize)

tmprs=rs
excllor2=F
if(excllor2){ # only plot cells with R^2>=rsqco
        tmprs[rs<rsqco]=NA
}

par(mfrow=c(2,3))
for(i in c(1:6)){
        xs(matrix(tmprs[,i],ncol=length(prsvals),nrow=length(latvals)),zlm=range(rs),maintxt=labels[i],units='r',lor2dots=F)
}
dev.off()

# for paper
png(paste(modeldir,'/Cross_Model_R_Conc_vs_Flux_SET_fig3a.png',sep=''),width=650,height=500*12/16,pointsize=12) # since colbar=T and colbarint=F gives heights=rep(c(12,2.0)
xs(matrix(tmprs[,1],ncol=length(prsvals),nrow=length(latvals)),zlm=max(abs(tmprs),na.rm=T)*c(-1,1),maintxt=expression('a. Southern Extratropical Flux Correlations'),lor2dots=F,colbar=F,color=blwtrd,legendarrows=F) # ,vtlines=c(-30,-10,0,20))
dev.off()
png(paste(modeldir,'/Cross_Model_R_Conc_vs_Flux_T_fig3b.png',sep=''),width=650,height=500,pointsize=12)
xs(matrix(tmprs[,2],ncol=length(prsvals),nrow=length(latvals)),zlm=max(abs(tmprs),na.rm=T)*c(-1,1),maintxt=expression('b. Tropical Flux Correlations'),units='correlation',lor2dots=F,colbarint=F,color=blwtrd,legendarrows=F) # ,vtlines=c(-30,-10,0,20))
dev.off()
png(paste(modeldir,'/Cross_Model_R_Conc_vs_Flux_NET_fig3c.png',sep=''),width=650,height=500*12/16,pointsize=12)
xs(matrix(tmprs[,3],ncol=length(prsvals),nrow=length(latvals)),zlm=max(abs(tmprs),na.rm=T)*c(-1,1),maintxt=expression('c. Northern Extratopical Flux Correlations'),units='r',lor2dots=F,colbar=F,color=blwtrd,legendarrows=F) # ,vtlines=c(-30,-10,0,20))
dev.off()
system(paste('montage ',modeldir,'/Cross_Model_R_Conc_vs_Flux_SET_fig3a.png ',modeldir,'/Cross_Model_R_Conc_vs_Flux_T_fig3b.png ',modeldir,'/Cross_Model_R_Conc_vs_Flux_NET_fig3c.png  -geometry +1+3 -gravity North ',modeldir,'/fig3_correlations.png',sep=''))

## EPS
width_cm <- 17.8/3
width_in <- width_cm / 2.54
sc <- 0.21*2/3 # Scaling factor for absolute line widths
new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)

file_name <- paste(modeldir,'/Cross_Model_R_Conc_vs_Flux_SET_fig3a.eps',sep='')
cairo_ps(filename = file_name,
         width = width_in,
         height = width_in*500/650*12/16,
         pointsize = new_pointsize)
xs(matrix(tmprs[,1],ncol=length(prsvals),nrow=length(latvals)),zlm=max(abs(tmprs),na.rm=T)*c(-1,1),maintxt=expression('a. Southern Extratropical Flux Correlations'),lor2dots=F,colbar=F,color=blwtrd,legendarrows=F,axlwd=0.3) # ,vtlines=c(-30,-10,0,20))
dev.off()
system(paste('cp ',modeldir,'/Cross_Model_R_Conc_vs_Flux_SET_fig3a.eps ',modeldir,'/fig3a.eps',sep=''))

file_name <- paste(modeldir,'/Cross_Model_R_Conc_vs_Flux_T_fig3b.eps',sep='')
cairo_ps(filename = file_name,
         width = width_in,
         height = width_in*500/650,
         pointsize = new_pointsize)
xs(matrix(tmprs[,2],ncol=length(prsvals),nrow=length(latvals)),zlm=max(abs(tmprs),na.rm=T)*c(-1,1),maintxt=expression('b. Tropical Flux Correlations'),units='correlation',lor2dots=F,colbarint=F,color=blwtrd,legendarrows=F,axlwd=0.3) # ,vtlines=c(-30,-10,0,20))
dev.off()
system(paste('cp ',modeldir,'/Cross_Model_R_Conc_vs_Flux_T_fig3b.eps ',modeldir,'/fig3b.eps',sep=''))

file_name <- paste(modeldir,'/Cross_Model_R_Conc_vs_Flux_NET_fig3c.eps',sep='')
cairo_ps(filename = file_name,
         width = width_in,
         height = width_in*500/650*12/16,
         pointsize = new_pointsize)
xs(matrix(tmprs[,3],ncol=length(prsvals),nrow=length(latvals)),zlm=max(abs(tmprs),na.rm=T)*c(-1,1),maintxt=expression('c. Northern Extratopical Flux Correlations'),units='r',lor2dots=F,colbar=F,color=blwtrd,legendarrows=F,axlwd=0.3) # ,vtlines=c(-30,-10,0,20))
dev.off()
system(paste('cp ',modeldir,'/Cross_Model_R_Conc_vs_Flux_NET_fig3c.eps ',modeldir,'/fig3c.eps',sep=''))


## look at T-ET corr
png(paste(modeldir,'/Cross_Model_R_Conc_vs_Flux_T-ET.png',sep=''),width=650,height=500,pointsize=12)
xs(matrix(tmprs[,5],ncol=length(prsvals),nrow=length(latvals)),zlm=max(abs(tmprs),na.rm=T)*c(-1,1),maintxt=expression('Tropical-Extratropical Flux Correlations'),units='correlation',lor2dots=F,colbarint=F,color=blwtrd,legendarrows=F)
dev.off()


# plot coefs
png(paste(modeldir,'/Cross_Model_Coefs_Conc_vs_Flux_byreg.png',sep=''),width=pngw*3,height=pngh*2,pointsize=psize)

tmpcoefs=coefs
if(excllor2){ 
	tmpcoefs[rsqds<rsqco]=NA
} 

par(mfrow=c(2,3))
for(i in c(1:6)){
	xs(matrix(tmpcoefs[,i],ncol=length(prsvals),nrow=length(latvals)),zlm=c(-1,1)*max(abs(coefs[,i])),maintxt=labels[i],color=blylrd,units='ppm/(PgC/yr)',lor2dots=F)
}
dev.off()

# plot ressd
png(paste(modeldir,'/Cross_Model_ResidualSDs_Conc_vs_Flux_byreg.png',sep=''),width=pngw*3,height=pngh*2,pointsize=psize)
tmpressd=ressd
if(excllor2){
        tmpressd[rsqds<rsqco]=NA
}

par(mfrow=c(2,3))
for(i in c(1:6)){
	xs(matrix(tmpressd[,i],ncol=length(prsvals),nrow=length(latvals)),maintxt=labels[i],lor2dots=F,zlm=range(c(0,as.vector(tmpressd)),na.rm=T))
}
dev.off()


# for paper
png(paste(modeldir,'/figS3_sdandeof1_ab.png',sep=''),width=650*2,height=500,pointsize=15)
xs(sdmat,maintxt=expression('a. Standard Deviations across Inversions'),zlm=range(c(0,as.vector(sdmat)),na.rm=T),colbarint=F,columns=2,units=expression(paste(CO[2],' (ppm)')),legendarrows=F)
i=1
xs(eof1,maintxt=expression('b. Empirical Orthogonal Function 1'),zlm=max(abs(eof1))*c(-1,1),color=blwtrd,colbarint=F,newlayout=F,units='',legendarrows=F) # ,zlm=c(-1,1)*max(abs(eof1))
dev.off()
system(paste('montage ',modeldir,'/figS3_sdandeof1_ab.png ',modeldir,'/figS3_sdandeof1_cd.png -geometry +1+2 -tile 1x2 ',modeldir,'/figS3_sdandeof1.png',sep=''))
