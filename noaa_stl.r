stlsite=function(site='mlo'){

	print(paste('Calculating STL fit for',site))

	infile=paste('LOCALDATA/co2_',site,'_surface-flask_1_ccgg_month.txt',sep='') # from e.g.: https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/txt/co2_mlo_surface-flask_1_ccgg_month.txt
	header=readLines(infile,n=1); hlines=as.numeric(tail(strsplit(header,' ')[[1]],1))
	indata=read.table(infile,skip=hlines,header=F,stringsAsFactors=F)
	colnames(indata)=c('site','year','month','value')
	plotdat=indata[,c('year','month','value')]; colnames(plotdat)=c('yr','mon','co2')
	plotdat$co2[plotdat$co2==-999.99]=NA; plotdat$co2[substr(indata$qcflag,1,1)!='.']=NA; plotdat$co2[substr(indata$qcflag,2,2)!='.']=NA
	plotdat=plotdat[plotdat$yr>1983,] # gets rid of all NAs, but just in case of new ones:
	plotdat$co2fill=plotdat$co2
	if(any(is.na(plotdat$co2fill))) plotdat$co2fill[is.na(plotdat$co2fill)]=approx(plotdat$decdate,plotdat$co2fill,plotdat$decdate[is.na(plotdat$co2fill)])$y # fills in any missing (bounded) values (does not extrapolate)
	plotdatetime=strptime(paste(plotdat$yr,plotdat$mon,'15'),format='%Y %m %d',tz='UTC')

	plotdat$date=paste(toupper(month.abb[plotdat$mon]),plotdat$yr,sep='')

	data.ts=ts(plotdat$co2fill,c(plotdatetime[1]$year+1900,plotdatetime[1]$mon+1),frequency=12)
	data.stl=stl(data.ts,s.window=5,t.window=121)
	co2dt=plotdat$co2fill-data.stl$time.series[,'trend']
	data.ts2=ts(co2dt,c(plotdatetime[1]$year+1900,plotdatetime[1]$mon+1),frequency=12)
	data.stl2=stl(data.ts2,s.window=5,t.window=25)

	bitmap(paste(site,'_co2_stl_%d.png',sep=''),width=10.5,height=7,res=100,type='png16',pointsize=10) # open device for creating png

	par(mfrow=c(1,1))
	par(oma=c(0,0,2,0))
	par(mar=c(3,3,2,1)+.1)
	par(mgp=c(2,1,0))

	plot(plotdatetime,plotdat$co2,type='n',ylab='ppm',main=toupper(site))
	points(plotdatetime,plotdat$co2fill,col='Dark Blue')
	points(plotdatetime,plotdat$co2,col='Red')

	par(mfrow=c(5,1))
	par(oma=c(5,0,4,0))
	par(mar=c(0,4,0,4))
	par(mgp=c(2.5,1,0))

	plot(plotdatetime,plotdat$co2,type='n',ylab='ppm',axes=F,cex.axis=1.5,cex.lab=1.5)
	box()
	axis(2,cex.axis=1.5)
	points(plotdatetime,plotdat$co2fill,col='Red')
	lines(plotdatetime,data.stl$time.series[,'trend']+data.stl2$time.series[,'trend']+data.stl2$time.series[,'seasonal'],col='Black')
	text(par('usr')[1],par('usr')[4],'a',adj=c(-1,2),cex=2)
	mtext(toupper(site),3,1,outer=T)

	plot(plotdatetime,data.stl$time.series[,'trend'],type='n',axes=F,ylab='',cex.axis=1.5,cex.lab=1.5)
	box()
	axis(4,cex.axis=1.5)
	mtext('ppm',4,2.5,cex=1.0)
	text(par('usr')[1],par('usr')[4],'b',adj=c(-1,2),cex=2)
	lines(plotdatetime,data.stl$time.series[,'trend'],col='Dark Green')

	plot(plotdatetime,data.stl2$time.series[,'trend'],type='n',ylab='ppm',axes=F,cex.axis=1.5,cex.lab=1.5,ylim=c(-1,1))
	box()
	axis(2,cex.axis=1.5)
	text(par('usr')[1],par('usr')[4],'c',adj=c(-1,2),cex=2)
	lines(plotdatetime,data.stl2$time.series[,'trend'],col='Dark Blue')

	plot(plotdatetime,data.stl2$time.series[,'seasonal'],type='n',axes=F,ylab='',cex.axis=1.5,cex.lab=1.5,ylim=c(-1.5,1.5))
	box()
	axis(4,cex.axis=1.5)
	mtext('ppm',4,2.5,cex=1.0)
	text(par('usr')[1],par('usr')[4],'d',adj=c(-1,2),cex=2)
	lines(plotdatetime,data.stl2$time.series[,'seasonal'],col='Green')

	plot(plotdatetime,data.stl2$time.series[,'remainder'],type='n',ylab='ppm',xlab='Year',cex.axis=1,cex.axis=1.5,cex.lab=1.5,ylim=c(-1,1))
	text(par('usr')[1],par('usr')[4],'e',adj=c(-1,2),cex=2)
	lines(plotdatetime,data.stl2$time.series[,'remainder'],col='Dark Blue')

	mtext('Year',1,3,outer=T,cex=1.2)

	par(mfrow=c(4,1))
	par(oma=c(5,0,4,0))
	par(mar=c(0,4,0,4))
	par(mgp=c(2.5,1,0))


	plot(plotdatetime,data.stl$time.series[,'trend'],type='n',axes=F,cex.axis=1.5,cex.lab=1.5,ylab=expression(paste(CO[2],' ',(ppm))))
	box()
	axis(2,cex.axis=1.5)
	text(par('usr')[1],par('usr')[4],'a',adj=c(-1,2),cex=2)
	lines(plotdatetime,data.stl$time.series[,'trend'],col='Red')
	mtext(toupper(site),3,1,outer=T)

	plot(plotdatetime,data.stl2$time.series[,'trend'],type='n',axes=F,cex.axis=1.5,cex.lab=1.5,ylim=c(-1,1),ylab='')
	box()
	axis(4,cex.axis=1.5)
	mtext(expression(paste(Delta,' ',CO[2],' ',(ppm))),4,2.5,cex=1.0)
	text(par('usr')[1],par('usr')[4],'b',adj=c(-1,2),cex=2)
	lines(plotdatetime,data.stl2$time.series[,'trend'],col='Dark Blue')

	plot(plotdatetime,data.stl2$time.series[,'seasonal'],type='n',axes=F,cex.axis=1.5,cex.lab=1.5,ylim=c(-1.5,1.5),ylab=expression(paste(Delta,' ',CO[2],' ',(ppm))))
	box()
	axis(2,cex.axis=1.5)
	text(par('usr')[1],par('usr')[4],'c',adj=c(-1,2),cex=2)
	lines(plotdatetime,data.stl2$time.series[,'seasonal'],col='Green')

	plot(plotdatetime,data.stl2$time.series[,'remainder'],type='n',cex.axis=1,cex.axis=1.5,cex.lab=1.5,ylim=c(-1,1),xlab='Year',ylab='',yaxt='n')
	axis(4,cex.axis=1.5)
	mtext(expression(paste(Delta,' ',CO[2],' ',(ppm))),4,2.5,cex=1.0)
	text(par('usr')[1],par('usr')[4],'d',adj=c(-1,2),cex=2)
	lines(plotdatetime,data.stl2$time.series[,'remainder'],col='Purple')

	dev.off()

	write('date co2 co2fill trend iav seasonal remainder total',paste(site,'_co2_stl_results.txt',sep=''))
	write(rbind(plotdat$date,round(rbind(plotdat$co2,plotdat$co2fill,data.stl$time.series[,'trend'],data.stl2$time.series[,'trend'],data.stl2$time.series[,'seasonal'],data.stl2$time.series[,'remainder'],data.stl$time.series[,'trend']+data.stl2$time.series[,'trend']+data.stl2$time.series[,'seasonal']),3)),paste(site,'_co2_stl_results.txt',sep=''),append=T,ncol=8)

} # end of stlsiteplot function
