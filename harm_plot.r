## program to plot the results of harmonic fits to binned data

harmp=function(species='CO2_NOAA',specname='CO2',units='pppm',surfco=900,coltop=300,detrend=T,dtrsta='mlo',xlim1=c(-70,90),ylim1=c(1000,0),zlimseas=NULL,zlimann=NULL,zlimphs=NULL,zlimcnt=NULL,zlimmon=NULL,inputdir='.',netdir='.',xtick=20,colbar=T,colbarint=F,pngw=650,pngh=650,psize=12,flag='',meanlegloc='bottomleft',amplegloc='bottomright',addthetacont=F,maintf=F,maintxt=NULL,legendarrows=F,truncextremes=F,type='img',out='png',model=NULL,tracer=NULL){ 

	# blue-yellow-red
	pal <- colorRampPalette(c('dark blue','blue','cyan','yellow','orange','red','red'))
	color2 <- pal(1000)
	colramp=colorRampPalette(color2)

	if(!is.null(model)){
		inputdir=paste(model,'/',tracer,sep='')
		netdir=inputdir ## this is just copying the files onto themselves presently
	}

	# read in harm_fit.r seasonal amplitude output 
	latvals=scan(paste(inputdir,'/Seasonal_Amplitude_',species,flag,'.txt',sep=''),nlines=1)
	prsvals=scan(paste(inputdir,'/Seasonal_Amplitude_',species,flag,'.txt',sep=''),nlines=1,skip=1)
	seasamp=read.table(paste(inputdir,'/Seasonal_Amplitude_',species,flag,'.txt',sep=''),skip=2)
	colseasamp=read.table(paste(inputdir,'/Seasonal_Amplitude_',species,flag,'_Column.txt',sep=''))
	surfseasamp=read.table(paste(inputdir,'/Seasonal_Amplitude_',species,flag,'_Surface.txt',sep=''))

	# flip pressure scale
	seasamp=seasamp[,ncol(seasamp):1]
	ylim1=1000-ylim1
	prsvals=rev(1000-prsvals)

	if(is.null(zlimseas)){ zlimseas=range(seasamp[,rev(prsvals)>(1000-ylim1[2])],na.rm=T) } 

	if(truncextremes){
		seasamp[seasamp<zlimseas[1]]=zlimseas[1]
		seasamp[seasamp>zlimseas[2]]=zlimseas[2]
	}

	if(addthetacont){
		thetaannmean=read.table(paste('Annual_Mean_THETA',flag,'.txt',sep=''),skip=2)
		thetaannmean=thetaannmean[,ncol(thetaannmean):1]
	}

	if(colbar){
		if(colbarint){
			#SCW implementation of cross section plot
			ylim2=c(ylim1[1],ylim1[1] + diff(ylim1)*1.30) # expand by 20% to accomodate legend
			yl1=ylim1[1] + diff(ylim1)*1.15;yl2=ylim1[1] + diff(ylim1)*1.25;
		} else { # external
			ylim2=ylim1
		}
		xl1=zlimseas[1];xl2=zlimseas[2]; XX=seq(xl1,xl2,length.out=length(color2));
		#XMAT is the matrix used for the image of the color bar
		XMAT=t(matrix(XX[2:length(XX)]-(XX[2]-XX[1])/2,nrow=1,byrow=T) )
	} else { # no color bar
		ylim2=ylim1
	}


	if(out=='png'){
		png(paste(netdir,'/Seasonal_Amplitude_',species,flag,'_xsect.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	} else if(out=='pdf'){
		pdf(paste(netdir,'/Seasonal_Amplitude_',species,flag,'_xsect.pdf',sep=''),width=pngw/650*7,height=pngh/650*7,pointsize=psize*0.75)
	}

	# plot main cross section image 

	if(maintf){
		par(mar=c(4, 4, 4, 0.5)+0.1)
		if(is.null(maintxt)){
			maintxt=paste('HIPPO/ORCAS/ATom Seasonal ',specname,' Amplitudes',sep='')
		} # else use input maintxt
	} else {
		par(mar=c(4, 4, 0, 0.5)+0.1)
		maintxt=''
		
	}

	xlb=expression(paste("Latitude (",degree, "N)")) 
	ylb="Pressure (mbar)"

	if(colbar&!colbarint){
		layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(12,2.5))
	}

	if(type=='img'){
		image(latvals,prsvals,as.matrix(seasamp),zlim=zlimseas,xlim=xlim1,ylim=ylim2,col=color2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F)
	}else if(type=='cnt'){
		plot(latvals,seq(prsvals[1],prsvals[2],length.out=length(latvals)),xlim=xlim1,ylim=ylim2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F,type='n',xaxs='i',yaxs='i')
		.filled.contour(latvals,prsvals,as.matrix(seasamp),levels = seq(zlimseas[1],zlimseas[2],length.out=20), col=colramp(19))
	}

	if(addthetacont){
		contour(latvals,prsvals,as.matrix(thetaannmean),add=T,col='gray35',levels=seq(270,370,10),vfont=c("serif", "bold"),lwd=3,labcex=0.9)
	}
	print(paste('range of amplitudes:',paste(range(seasamp,na.rm=T),collapse=', ')))
	print('5, 10, 50, 90, 95 % quantiles:')
	print(quantile(seasamp,c(0.05,0.1,0.5,0.9,0.95),na.rm=T))
	print(paste('lowest 5 of',sum(!is.na(seasamp))))
	print(seasamp[!is.na(seasamp)][order(seasamp[!is.na(seasamp)])][1:5])
	print(paste('highest 5 of',sum(!is.na(seasamp))))
	print(seasamp[!is.na(seasamp)][order(seasamp[!is.na(seasamp)],decreasing=T)][1:5]) 
	box()
	axis(1,at=seq(-100,100,xtick),labels=seq(-100,100,xtick),cex.axis=1.5)
	axis(2,at=1000-seq(1000,200,-200),labels=seq(1000,200,-200),cex.axis=1.5)

	if(colbar){
		if(colbarint){
			mtext(units,3,-8,cex=1.5)
			par(new=T); 
			#plot color bar
			image(x=XX,y=c(yl1,yl2),z=XMAT,zlim=zlimseas,xlim=c(zlimseas[1]-.1*diff(zlimseas),zlimseas[2]+.1*diff(zlimseas)),ylim=ylim2,axes=F,xlab="",ylab="",col=color2)
			lines(c(xl1,xl2,xl2,xl1,xl1),c(yl2,yl2,yl1,yl1,yl2))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=yl1,at=AA,cex.axis=1.5)
			axis(side=1,pos=yl1,at=AA,labels=F,tck=.075)
		} else { # colbarint=F (external)
			par(mar=c(3,5,1,1)) # impact of margins depends on height and width and psize of main plot
			par(mgp=c(1.7,1,0))
			squashhz=1.2
			squashvt=0.7
			image(x=XX,y=c(0,1),z=XMAT,zlim=zlimseas,xlim=c(zlimseas[1]-squashhz*diff(zlimseas),zlimseas[2]+squashhz*diff(zlimseas)),ylim=c(-squashvt,1+squashvt),axes=F,ylab="",col=color2,xlab=units,cex.lab=1.5)
			if(legendarrows){
				arrowwid=(xl2-xl1)*0.05
				polygon(c(xl2,xl2+arrowwid,xl2,xl2),c(0,0.5,1,0),col=tail(color2))
				polygon(c(xl1,xl1-arrowwid,xl1,xl1),c(0,0.5,1,0),col=head(color2,1))
				lines(c(xl2,xl2+arrowwid,xl2),c(0,0.5,1))
				lines(c(xl1,xl1-arrowwid,xl1),c(0,0.5,1))
			}
			lines(c(xl1,xl2,xl2,xl1,xl1),c(1,1,0,0,1))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=0,at=AA,cex.axis=1.5)
		} # int or ext
	} # if(colbar)

	dev.off()


	if(out=='png'){
		png(paste(netdir,'/Seasonal_Amplitude_',species,flag,'_vlat.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	} else if(out=='pdf'){
		pdf(paste(netdir,'/Seasonal_Amplitude_',species,flag,'_vlat.pdf',sep=''),width=pngw/650*7,height=pngh/650*7,pointsize=psize*0.75)
	}

	par(mar=c(5,5,4,2)+0.1)
	# plot column and surface
	ylm=range(c(colseasamp[,2],surfseasamp[,2]),na.rm=T)
	if(any(!is.na(colseasamp[,2]))){
		plot(latvals,colseasamp[,2],type='b',pch=21,col='blue',bg='white',xlab=xlb,ylab=paste(specname,' (',units,')',sep=''),main=maintxt,lwd=4,cex=2,cex.axis=2,cex.lab=2,cex.main=2,ylim=ylm,xlim=xlim1)
		points(latvals,surfseasamp[,2],type='b',pch=22,col='red',bg='white',lwd=4,cex=2)
		legend(amplegloc,c(paste('Column (>',coltop,')',sep=''),paste('Surface (>',surfco,')',sep='')),col=c('blue','red'),pch=c(21,22),pt.bg='white',lwd=4,pt.cex=2,cex=2)
	} else {
		plot(1:3,1:3,type='n',axes=F,xlab='',ylab='',main='No non-NA colseasamp')
	}

	dev.off()

	# Amplitude ratio plot
	if(out=='png'){
		png(paste(netdir,'/Seasonal_Amplitude_Ratio_',species,flag,'_vlat.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	} else if(out=='pdf'){
		pdf(paste(netdir,'/Seasonal_Amplitude_Ratio_',species,flag,'_vlat.pdf',sep=''),width=pngw/650*7,height=pngh/650*7,pointsize=psize*0.75)
	}

	par(mar=c(5,5,4,2)+0.1)
	# plot surface / column
	ylm=range(surfseasamp[,2]/colseasamp[,2],na.rm=T)
	if(any(!is.na(surfseasamp[,2]/colseasamp[,2]))){
		plot(latvals,surfseasamp[,2]/colseasamp[,2],type='b',pch=21,col='purple',bg='white',xlab=xlb,ylab=specname,main=maintxt,lwd=4,cex=2,cex.axis=2,cex.lab=2,cex.main=2,ylim=ylm,xlim=xlim1)
		legend(amplegloc,paste('Surface (>',surfco,') / Column (>',coltop,') Amplitude Ratio',sep=''),col=c('purple'),pch=c(21),pt.bg='white',lwd=4,pt.cex=2,cex=2)
	} else {
		plot(1:3,1:3,type='n',axes=F,xlab='',ylab='',main='No non-NA surfseasamp/colseasamp')
	}

	dev.off()


	# read in harm_fit.r annual mean output 
	latvals=scan(paste(inputdir,'/Annual_Mean_',species,flag,'.txt',sep=''),nlines=1)
	prsvals=scan(paste(inputdir,'/Annual_Mean_',species,flag,'.txt',sep=''),nlines=1,skip=1)
	annmean=read.table(paste(inputdir,'/Annual_Mean_',species,flag,'.txt',sep=''),skip=2)
	colannmean=read.table(paste(inputdir,'/Annual_Mean_',species,flag,'_Column.txt',sep=''))
	surfannmean=read.table(paste(inputdir,'/Annual_Mean_',species,flag,'_Surface.txt',sep=''))
	# flip pressure scale
	annmean=annmean[,ncol(annmean):1]
	prsvals=rev(1000-prsvals)

	if(is.null(zlimann)){ zlimann=range(annmean[,rev(prsvals)>(1000-ylim1[2])],na.rm=T) } 

	if(truncextremes){
		annmean[annmean<zlimann[1]]=zlimann[1]
		annmean[annmean>zlimann[2]]=zlimann[2]
	}


	xl1=zlimann[1];xl2=zlimann[2]; XX=seq(xl1,xl2,length.out=length(color2));
	#XMAT is the matrix used for the image of the color bar
	XMAT=t(matrix(XX[2:length(XX)]-(XX[2]-XX[1])/2,nrow=1,byrow=T) )

	# plot main cross section image 

	if(out=='png'){
		png(paste(netdir,'/Annual_Mean_',species,flag,'_xsect.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	} else if(out=='pdf'){
		pdf(paste(netdir,'/Annual_Mean_',species,flag,'_xsect.pdf',sep=''),width=pngw/650*7,height=pngh/650*7,pointsize=psize*0.75)
	}

	if(maintf){
		par(mar=c(4, 4, 4, 0.5)+0.1)
		if(is.null(maintxt)){
			maintxt=paste('HIPPO/ORCAS/ATom Annual ',specname,' Means',sep='')
		} # else use input maintxt
	} else {
		par(mar=c(4, 4, 0, 0.5)+0.1)
		maintxt=''

	}

	xlb=expression(paste("Latitude (",degree, "N)"))
	ylb="Pressure (mbar)"

	if(colbar&!colbarint){
		layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(12,2.5))
	}

	if(type=='img'){
		image(latvals,prsvals,as.matrix(annmean),zlim=zlimann,xlim=xlim1,ylim=ylim2,col=color2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F)
	}else if(type=='cnt'){
		plot(latvals,seq(prsvals[1],prsvals[2],length.out=length(latvals)),xlim=xlim1,ylim=ylim2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F,type='n',xaxs='i',yaxs='i')
		.filled.contour(latvals,prsvals,as.matrix(annmean),levels = seq(zlimann[1],zlimann[2],length.out=20), col=colramp(19))
	}

	if(addthetacont){
		contour(latvals,prsvals,as.matrix(thetaannmean),add=T,col='gray35',levels=seq(270,370,10),vfont=c("serif", "bold"),lwd=3,labcex=0.9)
	}
	print(paste('range of annual means:',paste(range(annmean,na.rm=T),collapse=', ')))
	print('5, 10, 50, 90, 95 % quantiles:')
	print(quantile(annmean,c(0.05,0.1,0.5,0.9,0.95),na.rm=T))
	print(paste('lowest 5 of',sum(!is.na(annmean))))
	print(annmean[!is.na(annmean)][order(annmean[!is.na(annmean)])][1:5])
	print(paste('highest 5 of',sum(!is.na(annmean))))
	print(annmean[!is.na(annmean)][order(annmean[!is.na(annmean)],decreasing=T)][1:5]) 
	box()
	axis(1,at=seq(-100,100,xtick),labels=seq(-100,100,xtick),cex.axis=1.5)
	axis(2,at=1000-seq(1000,200,-200),labels=seq(1000,200,-200),cex.axis=1.5)

	if(colbar){
		if(colbarint){
			mtext(units,3,-8,cex=1.5)
			par(new=T);
			#plot color bar
			image(x=XX,y=c(yl1,yl2),z=XMAT,zlim=zlimann,xlim=c(zlimann[1]-.1*diff(zlimann),zlimann[2]+.1*diff(zlimann)),ylim=ylim2,axes=F,xlab="",ylab="",col=color2)
			lines(c(xl1,xl2,xl2,xl1,xl1),c(yl2,yl2,yl1,yl1,yl2))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=yl1,at=AA,cex.axis=1.5)
			axis(side=1,pos=yl1,at=AA,labels=F,tck=.075)
		} else { # colbarint=F (external)
			par(mar=c(3,5,1,1)) # impact of margins depends on height and width and psize of main plot
			par(mgp=c(1.7,1,0))
			squashhz=1.2
			squashvt=0.7
			image(x=XX,y=c(0,1),z=XMAT,zlim=zlimann,xlim=c(zlimann[1]-squashhz*diff(zlimann),zlimann[2]+squashhz*diff(zlimann)),ylim=c(-squashvt,1+squashvt),axes=F,ylab="",col=color2,xlab=units,cex.lab=1.5)
			if(legendarrows){
				arrowwid=(xl2-xl1)*0.05
				polygon(c(xl2,xl2+arrowwid,xl2,xl2),c(0,0.5,1,0),col=tail(color2))
				polygon(c(xl1,xl1-arrowwid,xl1,xl1),c(0,0.5,1,0),col=head(color2,1))
				lines(c(xl2,xl2+arrowwid,xl2),c(0,0.5,1))
				lines(c(xl1,xl1-arrowwid,xl1),c(0,0.5,1))
			}
			lines(c(xl1,xl2,xl2,xl1,xl1),c(1,1,0,0,1))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=0,at=AA,cex.axis=1.5)
		} # int or ext
	} # if(colbar)

	dev.off()

	if(out=='png'){
		png(paste(netdir,'/Annual_Mean_',species,flag,'_vlat.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	} else if(out=='pdf'){
		pdf(paste(netdir,'/Annual_Mean_',species,flag,'_vlat.pdf',sep=''),width=pngw/650*7,height=pngh/650*7,pointsize=psize*0.75)
	}

	par(mar=c(5,5,4,2)+0.1)
	# plot column and surface
	ylm=range(c(colannmean[,2],surfannmean[,2]),na.rm=T)
	if(detrend){ 
		ylb=paste(specname,' (',units,') - ',toupper(dtrsta),' trend',sep='')
	} else {
		ylb=paste(specname,' (',units,')',sep='')
	}

	if(any(!is.na(colannmean[,2]))){
		plot(latvals,colannmean[,2],type='b',pch=21,col='blue',bg='white',xlab=xlb,ylab=ylb,main=maintxt,lwd=4,cex=2,cex.axis=2,cex.lab=2,cex.main=2,ylim=ylm,xlim=xlim1)
		points(latvals,surfannmean[,2],type='b',pch=22,col='red',bg='white',lwd=4,cex=2)
		legend(meanlegloc,c(paste('Column (>',coltop,' mbar)',sep=''),paste('Surface (>',surfco,' mbar)',sep='')),col=c('blue','red'),pch=c(21,22),pt.bg='white',lwd=4,pt.cex=2,cex=2)
	} else {
		plot(1:3,1:3,type='n',axes=F,xlab='',ylab='',main='No non-NA colseasamp')
	}

	dev.off()

	# read in harm_fit.r phase output 
	latvals=scan(paste(inputdir,'/Seasonal_Phase_',species,flag,'.txt',sep=''),nlines=1)
	prsvals=scan(paste(inputdir,'/Seasonal_Phase_',species,flag,'.txt',sep=''),nlines=1,skip=1)
	phase=read.table(paste(inputdir,'/Seasonal_Phase_',species,flag,'.txt',sep=''),skip=2)
	colphase=read.table(paste(inputdir,'/Seasonal_Phase_',species,flag,'_Column.txt',sep=''))
	surfphase=read.table(paste(inputdir,'/Seasonal_Phase_',species,flag,'_Surface.txt',sep=''))
	# flip pressure scale
	phase=phase[,ncol(phase):1]
	prsvals=rev(1000-prsvals)

	if(is.null(zlimphs)){ zlimphs=range(phase[,rev(prsvals)>(1000-ylim1[2])],na.rm=T) } 

	xl1=zlimphs[1];xl2=zlimphs[2]; XX=seq(xl1,xl2,length.out=length(color2));

	#XMAT is the matrix used for the image of the color bar
	XMAT=t(matrix(XX[2:length(XX)]-(XX[2]-XX[1])/2,nrow=1,byrow=T) )

	# plot main cross section image 

	if(out=='png'){
		png(paste(netdir,'/Seasonal_Phase_',species,flag,'_xsect.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	} else if(out=='pdf'){
		pdf(paste(netdir,'/Seasonal_Phase_',species,flag,'_xsect.pdf',sep=''),width=pngw/650*7,height=pngh/650*7,pointsize=psize*0.75)
	}

	if(maintf){
		par(mar=c(4, 4, 4, 0.5)+0.1)
		if(is.null(maintxt)){
			maintxt=paste('HIPPO/ORCAS/ATom Seasonal ',specname,' Phase',sep='')
		} # else use input maintxt
	} else {
		par(mar=c(4, 4, 0, 0.5)+0.1)
		maintxt=''

	}

	xlb=expression(paste("Latitude (",degree, "N)"))
	ylb="Pressure (mbar)"

	if(colbar&!colbarint){
		layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(12,2.5))
	}

	if(type=='img'){
		image(latvals,prsvals,as.matrix(phase),zlim=zlimphs,xlim=xlim1,ylim=ylim2,col=color2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F)
	}else if(type=='cnt'){
		plot(latvals,seq(prsvals[1],prsvals[2],length.out=length(latvals)),xlim=xlim1,ylim=ylim2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F,type='n',xaxs='i',yaxs='i')
		.filled.contour(latvals,prsvals,as.matrix(phase),levels = seq(zlimphs[1],zlimphs[2],length.out=20), col=colramp(19))
	}

	print(paste('range of seasonal phase:',paste(range(phase,na.rm=T),collapse=', ')))
	print('5, 10, 50, 90, 95 % quantiles:')
	print(quantile(phase,c(0.05,0.1,0.5,0.9,0.95),na.rm=T))
	print(paste('lowest 5 of',sum(!is.na(phase))))
	print(phase[!is.na(phase)][order(phase[!is.na(phase)])][1:5])
	print(paste('highest 5 of',sum(!is.na(phase))))
	print(phase[!is.na(phase)][order(phase[!is.na(phase)],decreasing=T)][1:5]) 
	box()
	axis(1,at=seq(-100,100,xtick),labels=seq(-100,100,xtick),cex.axis=1.5)
	axis(2,at=1000-seq(1000,200,-200),labels=seq(1000,200,-200),cex.axis=1.5)

	if(colbar){
		if(colbarint){
			mtext('Day of Minimum',3,-8,cex=1.5)
			par(new=T);
			#plot color bar
			image(x=XX,y=c(yl1,yl2),z=XMAT,zlim=zlimphs,xlim=c(zlimphs[1]-.1*diff(zlimphs),zlimphs[2]+.1*diff(zlimphs)),ylim=ylim2,axes=F,xlab="",ylab="",col=color2)
			lines(c(xl1,xl2,xl2,xl1,xl1),c(yl2,yl2,yl1,yl1,yl2))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=yl1,at=AA,cex.axis=1.5)
			axis(side=1,pos=yl1,at=AA,labels=F,tck=.075)
		} else { # colbarint=F (external)
			par(mar=c(3,5,1,1)) # impact of margins depends on height and width and psize of main plot
			par(mgp=c(1.7,1,0))
			squashhz=1.2
			squashvt=0.7
			image(x=XX,y=c(0,1),z=XMAT,zlim=zlimphs,xlim=c(zlimphs[1]-squashhz*diff(zlimphs),zlimphs[2]+squashhz*diff(zlimphs)),ylim=c(-squashvt,1+squashvt),axes=F,ylab="",col=color2,xlab='Day of Minimum',cex.lab=1.5)
			lines(c(xl1,xl2,xl2,xl1,xl1),c(1,1,0,0,1))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=0,at=AA,cex.axis=1.5)
		} # int or ext
	} # if(colbar)

	dev.off()

	# read in harm_fit.r count 
	latvals=scan(paste(inputdir,'/Count_',species,flag,'.txt',sep=''),nlines=1)
	prsvals=scan(paste(inputdir,'/Count_',species,flag,'.txt',sep=''),nlines=1,skip=1)
	count=read.table(paste(inputdir,'/Count_',species,flag,'.txt',sep=''),skip=2)
	# flip pressure scale
	count=count[,ncol(count):1]
	prsvals=rev(1000-prsvals)

	if(is.null(zlimcnt)){ zlimcnt=range(count[,rev(prsvals)>(1000-ylim1[2])],na.rm=T) } 

	xl1=zlimcnt[1];xl2=zlimcnt[2]; XX=seq(xl1,xl2,length.out=length(color2));

	#XMAT is the matrix used for the image of the color bar
	XMAT=t(matrix(XX[2:length(XX)]-(XX[2]-XX[1])/2,nrow=1,byrow=T) )

	# plot main cross section image 

	if(out=='png'){
		png(paste(netdir,'/Count_',species,flag,'_xsect.png',sep=''),width=pngw,height=pngh,pointsize=psize)
	} else if(out=='pdf'){
		pdf(paste(netdir,'/Count_',species,flag,'_xsect.pdf',sep=''),width=pngw/650*7,height=pngh/650*7,pointsize=psize*0.75)
	}

	if(maintf){
		par(mar=c(4, 4, 4, 0.5)+0.1)
		if(is.null(maintxt)){
			maintxt=paste('HIPPO/ORCAS/ATom ',specname,' Flight Count',sep='')
		} # else use input maintxt
	} else {
		par(mar=c(4, 4, 0, 0.5)+0.1)
		maintxt=''

	}

	xlb=expression(paste("Latitude (",degree, "N)"))
	ylb="Pressure (mbar)"

	if(colbar&!colbarint){
		layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(12,2.5))
	}

	if(type=='img'){
		image(latvals,prsvals,as.matrix(count),zlim=zlimcnt,xlim=xlim1,ylim=ylim2,col=color2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F)
	}else if(type=='cnt'){
		plot(latvals,seq(prsvals[1],prsvals[2],length.out=length(latvals)),xlim=xlim1,ylim=ylim2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F,type='n',xaxs='i',yaxs='i')
		.filled.contour(latvals,prsvals,as.matrix(count),levels = seq(zlimcnt[1],zlimcnt[2],length.out=20), col=colramp(19))
	}

	print(paste('range of seasonal count:',paste(range(count,na.rm=T),collapse=', ')))
	print('5, 10, 50, 90, 95 % quantiles:')
	print(quantile(count,c(0.05,0.1,0.5,0.9,0.95),na.rm=T))
	print(paste('lowest 5 of',sum(!is.na(count))))
	print(count[!is.na(count)][order(count[!is.na(count)])][1:5])
	print(paste('highest 5 of',sum(!is.na(count))))
	print(count[!is.na(count)][order(count[!is.na(count)],decreasing=T)][1:5]) 
	box()
	axis(1,at=seq(-100,100,xtick),labels=seq(-100,100,xtick),cex.axis=1.5)
	axis(2,at=1000-seq(1000,200,-200),labels=seq(1000,200,-200),cex.axis=1.5)

	if(colbar){
		if(colbarint){
			mtext('Flight Count',3,-8,cex=1.5)
			par(new=T);
			#plot color bar
			image(x=XX,y=c(yl1,yl2),z=XMAT,zlim=zlimcnt,xlim=c(zlimcnt[1]-.1*diff(zlimcnt),zlimcnt[2]+.1*diff(zlimcnt)),ylim=ylim2,axes=F,xlab="",ylab="",col=color2)
			lines(c(xl1,xl2,xl2,xl1,xl1),c(yl2,yl2,yl1,yl1,yl2))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=yl1,at=AA,cex.axis=1.5)
			axis(side=1,pos=yl1,at=AA,labels=F,tck=.075)
		} else { # colbarint=F (external)
			par(mar=c(3,5,1,1)) # impact of margins depends on height and width and psize of main plot
			par(mgp=c(1.7,1,0))
			squashhz=1.2
			squashvt=0.7
			image(x=XX,y=c(0,1),z=XMAT,zlim=zlimcnt,xlim=c(zlimcnt[1]-squashhz*diff(zlimcnt),zlimcnt[2]+squashhz*diff(zlimcnt)),ylim=c(-squashvt,1+squashvt),axes=F,ylab="",col=color2,xlab='Flight Count',cex.lab=1.5)
			lines(c(xl1,xl2,xl2,xl1,xl1),c(1,1,0,0,1))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=0,at=AA,cex.axis=1.5)
		} # int or ext
	} # if(colbar)

	dev.off()

	## read in coefficients and plot monthly mid-point snapshots
	coefs=read.table(paste(inputdir,'/Fit_Coefficients_',species,flag,'.txt',sep=''),header=T)
	nharm=(ncol(coefs)-3)/2 # set in harm_fit, this is just determining from dims
	mmid=seq(0.5,11.5)/12*2*pi #  ~monthly mid-points
	harm=cbind(rep(1,length(mmid)),cos(mmid),sin(mmid),cos(2*mmid),sin(2*mmid),cos(3*mmid),sin(3*mmid),cos(4*mmid),sin(4*mmid))[,1:(nharm*2+1)]

	mmidvals=array(NA,dim=c(length(latvals),length(prsvals),12))
	mmidanom=array(NA,dim=c(length(latvals),length(prsvals),12))
	print(latvals)
	print(prsvals)
	for(i in c(1:length(latvals))){
	for(j in c(1:length(prsvals))){
		if(any(!is.na(coefs[coefs$lat==i&coefs$prs==j,3:ncol(coefs)]))){
			mmidvals[i,j,]=drop(as.matrix(coefs[coefs$lat==i&coefs$prs==j,3:ncol(coefs)])%*%t(harm))
		}
	}
	}
	for(j in c(1:length(prsvals))){
			mmidanom[,j,]=sweep(mmidvals[,j,],1,colannmean[,2],"-")
	}

	if(out=='png'){
		png(paste(netdir,'/Monthly_Slices_',species,flag,'.png',sep=''),width=1800,height=1200,pointsize=24)
	} else if(out=='pdf'){
		pdf(paste(netdir,'/Monthly_Slices_',species,flag,'.pdf',sep=''),width=1800/650*7,height=1200/650*7,pointsize=psize*0.75)
	}

	par(mar=c(3.7, 4.2, 2, 0)+0.1)
	par(oma=c(1,1,1,1))

	layout(matrix(c(1:12,13,13,13,13),4,4,byrow=T),heights=c(6,6,6,3))

	# plot main cross section image
	print(quantile(mmidvals,probs=c(.95,.99),na.rm=T))
	if(is.null(zlimmon)){ zlimmon=range(mmidvals[,prsvals>(1000-ylim1[2]),],na.rm=T) } 
	for(i in (1:12)){
		print(quantile(mmidvals[,length(prsvals):1,i],probs=c(.95,.99),na.rm=T))
		maintxt=month.abb[i]
		if(is.element(i,c(1,5,9))){ ylb='Pressure (mbar)' } else { ylb='' }
		if(i>8){ xlb=expression(paste("Latitude (",degree, "N)")) } else { xlb='' }
		tmp=as.matrix(mmidvals[,length(prsvals):1,i])
		tmp[tmp>zlimmon[2]]=zlimmon[2]
		tmp[tmp<zlimmon[1]]=zlimmon[1]
		image(latvals,prsvals,tmp,zlim=zlimmon,xlim=xlim1,ylim=ylim1,col=color2,xlab=xlb,ylab=ylb,cex.axis=1.5,cex.lab=1.5,main=maintxt,cex.main=1.5,axes=F)
		box()
		axis(1,at=seq(-100,100,xtick),labels=seq(-100,100,xtick),cex.axis=1.5)
		axis(2,at=1000-seq(1000,200,-200),labels=seq(1000,200,-200),cex.axis=1.5)
	}

	# now add color bar
	zlm=zlimmon
	xl1=zlm[1];xl2=zlm[2]; XX=seq(xl1,xl2,length.out=length(color2)); # XX is edges of cells, XMAT is center
	XMAT=t(matrix(XX[2:length(XX)]-(XX[2]-XX[1])/2,nrow=1,byrow=T) ) #XMAT is the matrix used for the image of the color bar
	par(mar=c(3,1,1,1)) # impact of margins depends on height and width and psize of main plot
	par(mgp=c(1.7,1,0))
	squashhz=1.2
	squashvt=0.7
	image(x=XX,y=c(0,1),z=XMAT,zlim=zlm,xlim=c(zlm[1]-squashhz*diff(zlm),zlm[2]+squashhz*diff(zlm)),ylim=c(-squashvt,1+squashvt),axes=F,ylab="",col=color2,xlab=units,cex.lab=1.5)
	lines(c(xl1,xl2,xl2,xl1,xl1),c(1,1,0,0,1))
	breaks.axis=NULL # can specify color bar tick locations
	if(is.null(breaks.axis)){ aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)] }else{AA=breaks.axis}
	axis(side=1,pos=0,at=AA,cex.axis=1.5)

	dev.off()


	print('making index html')

	write('<HTML><TITLE> HIPPO/ORCAS/ATom HARMONIC FIT PAGE </TITLE><BODY TEXT="#000000" BGCOLOR="#FFFFFF">','index.html')
	dstamp<-system('date',intern=T)
	write('<H2>HIPPO/ORCAS/ATom FITS<BR>','index.html',append=T)
	write(paste('PROCESSED ',dstamp,'</H2>',sep=''),'index.html',append=T)
	write(paste('<A HREF="Seasonal_Amplitude_',species,flag,'_xsect.png"><img src="Seasonal_Amplitude_',species,flag,'_xsect.png" title="Seasonal Amplitudes" style="border: 0 px solid ; width: 205px; height: 120px;"></A></td>',sep=''),'index.html',append=T)
	write(paste('<A HREF="Seasonal_Amplitude_',species,flag,'_vlat.png"><img src="Seasonal_Amplitude_',species,flag,'_vlat.png" title="Column Seasonal Amplitudes" style="border: 0 px solid ; width: 205px; height: 120px;"></A></td>',sep=''),'index.html',append=T)
	write(paste('<A HREF="Annual_Mean_',species,flag,'_xsect.png"><img src="Annual_Mean_',species,flag,'_xsect.png" title="Annual Means" style="border: 0 px solid ; width: 205px; height: 120px;"></A></td>',sep=''),'index.html',append=T)
	write(paste('<A HREF="Annual_Mean_',species,flag,'_vlat.png"><img src="Annual_Mean_',species,flag,'_vlat.png" title="Annual Means" style="border: 0 px solid ; width: 205px; height: 120px;"></A></td>',sep=''),'index.html',append=T)
	write('<BR>','index.html',append=T)

	write('</BODY></HTML>','index.html',append=T)

	system(paste('mv index.html ',netdir,sep=''))
	system(paste('chmod 664 ',netdir,'/index.html',sep=''))

} # end of harmp function
