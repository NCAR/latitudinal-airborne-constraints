## define xsect plotting function
xs=function(zmat=allexparr[,,which(experiments=='IS')],zlm=NULL,maintxt='',xlb=expression(paste("Latitude (",degree, "N)")),ylb="Pressure (hPa)",units='',colbar=T,colbarint=T,legendarrows=T,color=usecol,lor2dots=F,columns=1,newlayout=T,vtlines=NULL,txtcex=1.5,axlwd=1){

	par(mar=c(4,5,4,1))
	par(mgp=c(3,1,0))

	if(is.null(zlm)) zlm=range(zmat)

	if(colbar){
		if(colbarint){
			ylim2=c(ylim1[1],ylim1[1] + diff(ylim1)*1.30) # expand by 20% to accomodate legend
			yl1=ylim1[1] + diff(ylim1)*1.15;yl2=ylim1[1] + diff(ylim1)*1.25;
        	} else { # external
                	ylim2=ylim1
        	}
		xl1=zlm[1];xl2=zlm[2]; XX=seq(xl1,xl2,length.out=length(color));
		XMAT=t(matrix(XX[2:length(XX)]-(XX[2]-XX[1])/2,nrow=1,byrow=T) )
	} else {
		ylim2=ylim1
	}

	if(colbar&!colbarint&newlayout){
        	layout(matrix(c(1:(2*columns)), nrow=2, ncol=columns), heights=rep(c(12,4.0),columns))
	}

	image(latvals,prsvals,zmat,zlim=zlm,xlim=xlim1,ylim=ylim2,col=color,xlab=xlb,ylab=ylb,cex.axis=txtcex,cex.lab=txtcex,main=maintxt,cex.main=txtcex,axes=F)
	box(lwd=axlwd)
	axis(1,at=seq(-100,100,xtick),labels=seq(-100,100,xtick),cex.axis=txtcex,lwd=axlwd)
	axis(2,at=seq(1000,200,-200),cex.axis=txtcex,lwd=axlwd)
	if(!is.null(vtlines)) abline(v=vtlines,lty='dashed',lwd=2,col='white')
	if(lor2dots) points(alllat[rsqds[,i]<rsqco],allprs[rsqds[,i]<rsqco],pch=16,cex=2,col='white')
	if(colbar){

		if(colbarint){

			mtext(units,3,-12,cex=txtcex)
			par(new=T);
			image(x=XX,y=c(yl2,yl1),z=XMAT,zlim=zlm,xlim=c(zlm[1]-.1*diff(zlm),zlm[2]+.1*diff(zlm)),ylim=ylim2,axes=F,xlab="",ylab="",col=color)
			lines(c(xl1,xl2,xl2,xl1,xl1),c(yl2,yl2,yl1,yl1,yl2))
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=yl1,at=AA,cex.axis=txtcex)
			axis(side=1,pos=yl1,at=AA,labels=F,tck=.075)

		} else { # colbarint=F (external)

			par(mar=c(3,5,0.1,1)) # impact of margins depends on height and width and psize of main plot
			par(mgp=c(1.3,1,0))
			squashhz=0.7
			squashvt=1.1
			image(x=XX,y=c(0,1),z=XMAT,zlim=zlm,xlim=c(zlm[1]-squashhz*diff(zlm),zlm[2]+squashhz*diff(zlm)),ylim=c(-squashvt,1+squashvt),axes=F,ylab="",col=color,xlab=units,cex.lab=txtcex)
			if(legendarrows){
				arrowwid=(xl2-xl1)*0.05
				polygon(c(xl2,xl2+arrowwid,xl2,xl2),c(0,0.5,1,0),col=tail(color,1),border=NA)
				polygon(c(xl1,xl1-arrowwid,xl1,xl1),c(0,0.5,1,0),col=head(color,1),border=NA)
				lines(c(xl2,xl2+arrowwid,xl2),c(0,0.5,1),lwd=axlwd)
				lines(c(xl1,xl1-arrowwid,xl1),c(0,0.5,1),lwd=axlwd)
				lines(c(xl1,xl2),c(1,1),lwd=axlwd)
				lines(c(xl1,xl2),c(0,0),lwd=axlwd)
			} else {
				lines(c(xl1,xl2,xl2,xl1,xl1),c(1,1,0,0,1),lwd=axlwd)
			}
			aa=pretty(XX,8);AA=aa[aa>=min(XX)&aa<=max(XX)]
			axis(side=1,pos=0,at=AA,cex.axis=txtcex,lwd=axlwd)

		} # int or ext

	}

}
