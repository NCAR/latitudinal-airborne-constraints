## program to plot ECFCs and return results

ecplotmlr<-function(xvar1=seq(1:6),xvar2=NULL,xvar3=NULL,yvar=c(1.5,3,2,5,4,6.5),obsval1=rep(5,11),obsval2=NULL,obsval3=NULL,filenames=c('ecplotout.png'),nsamp2d=1000,nsampmlr=1000,modnames=c('CT','CAMS','CSU','CT','CAMS','CSU'),expnames=c('IS','IS','IS','LNLG','LNLG','LNLG'),plotsel=c(1),fitexpmeans=F,inclmodleg=T,inclexpleg=T,xlb=expression('Concentration Metric (ppm)'),ylb=expression(paste('Total Flux (PgC ',yr^-1,')')),ecmain=expression('Emergent Constraint'),xlm=NULL,ylm=NULL,r2loc=NULL){
# example call
#ecplotmlr(xvar1=annnet[,i],xvar2-anntrop[,i],xvar3=NULL,yvar=regfluxes[,i],obsval1=obsannnet[,i],obsval2=obsanntrop[,i],obsval3=NULL,nsamp2d=100000,modnames=modnames,expnames=expnames,filename='NET_Flux_vs_NET-T_Conc_latdiv_OLS_Fig5c.png',plotsel=c(3),fitexpmeans=F,inclmodleg=F,inclexpleg=F,ecmain=expression('c. Northern Extratropics'),xlb=expression(paste('ATom Concentration Prediction (PgC ',yr^-     1,')')),ylb='',xlm=range(regfluxes[,i]),ylm=range(regfluxes[,i]) ) # ,ylb=expression(paste('v10 MIP Total Flux (PgC ',yr^-1,')')) )


	# define variables for primary flux estimate and obs uncertainty
	variables=names(obsval1) # c('CO2_OP','CO2_AO2_m_CO2_NOAA','CO2_QCLS_m_CO2_NOAA','CO2_MED_m_CO2_NOAA','CO2_PFP_m_CO2_NOAA')
	varsel='CO2_OP' # has to be in variables
	sdsel=c('CO2_AO2_m_CO2_NOAA','CO2_QCLS_m_CO2_NOAA','CO2_MED_m_CO2_NOAA','CO2_PFP_m_CO2_NOAA')

	library('MASS')
	library('vcov')

	# calc experiment means using only models with all experiments present
	expcount=aggregate(!is.na(xvar1),by=list(mod=modnames),sum)
	sel<-expcount$x[match(modnames,expcount$mod)]==length(unique(expnames))
	expmeanconc1=aggregate(xvar1[sel],by=list(exp=expnames[sel]),mean)$x
	expsdconc1=aggregate(xvar1[sel],by=list(exp=expnames[sel]),sd)$x
	if(!is.null(xvar2)){
		expmeanconc2=aggregate(xvar2[sel],by=list(exp=expnames[sel]),mean)$x
	        expsdconc2=aggregate(xvar2[sel],by=list(exp=expnames[sel]),sd)$x
		if(!is.null(xvar3)){ # can't have 3 without 2
			expmeanconc2=aggregate(xvar2[sel],by=list(exp=expnames[sel]),mean)$x
		        expsdconc2=aggregate(xvar2[sel],by=list(exp=expnames[sel]),sd)$x
		}
	}
	expmeanflux=aggregate(yvar[sel],by=list(exp=expnames[sel]),mean)$x
	expsdflux=aggregate(yvar[sel],by=list(exp=expnames[sel]),sd)$x

	# calculate a MLR fit to determine scaling on each region and covariation of parameters
	if(fitexpmeans){
		if(!is.null(xvar3)){
			mlrfit=lm(expmeanflux ~ expmeanconc1 + expmeanconc2 + expmeanconc3)
			print(summary(mlrfit)$coef)
		} else if(!is.null(xvar2)){
			mlrfit=lm(expmeanflux ~ expmeanconc1 + expmeanconc2)
		} else {
			mlrfit=lm(expmeanflux ~ expmeanconc1)
		}
	} else {
		if(!is.null(xvar3)){
			mlrfit=lm(yvar ~ xvar1 + xvar2 + xvar3)
			print(summary(mlrfit)$coef)
		} else if(!is.null(xvar2)){
			mlrfit=lm(yvar ~ xvar1 + xvar2)
		} else {
			mlrfit=lm(yvar ~ xvar1)
		}
	}

	covar=vcov(mlrfit)

	# for 2-D representation, calculate xvar, obsval, and expmeanconc
	if(!is.null(xvar3)){
		xvar2d=mlrfit$coef[1] + mlrfit$coef[2] * xvar1 + mlrfit$coef[3] * xvar2 + mlrfit$coef[4] * xvar3 
		obsval2d=mlrfit$coef[1] + mlrfit$coef[2] * obsval1[variables==varsel] + mlrfit$coef[3] * obsval2[variables==varsel] + mlrfit$coef[4] * obsval3[variables==varsel] 
	} else if(!is.null(xvar2)){
		xvar2d=mlrfit$coef[1] + mlrfit$coef[2] * xvar1 + mlrfit$coef[3] * xvar2
		obsval2d=mlrfit$coef[1] + mlrfit$coef[2] * obsval1[variables==varsel] + mlrfit$coef[3] * obsval2[variables==varsel]
       	} else {
		xvar2d=mlrfit$coef[1] + mlrfit$coef[2] * xvar1
		obsval2d=mlrfit$coef[1] + mlrfit$coef[2] * obsval1[variables==varsel]
       	}
	expmeanconc=aggregate(xvar2d[sel],by=list(exp=expnames[sel]),mean)$x
	expsdconc=aggregate(xvar2d[sel],by=list(exp=expnames[sel]),sd)$x

	# calculate obssd using mlrfit coef
	obsdifmet=0 # equiv of NOAA-NOAA
	for(difvar in sdsel){ # for sensor-NOAA
		if(!is.null(xvar3)){
        		obsdifmet=c(obsdifmet,mlrfit$coef[2] * obsval1[variables==difvar] + mlrfit$coef[3] * obsval2[variables==difvar] + mlrfit$coef[4] * obsval3[variables==difvar])
	        } else if(!is.null(xvar2)){
        		obsdifmet=c(obsdifmet,mlrfit$coef[2] * obsval1[variables==difvar] + mlrfit$coef[3] * obsval2[variables==difvar])
		} else {
        		obsdifmet=c(obsdifmet,mlrfit$coef[2] * obsval1[variables==difvar])
		}
	}
	obssd=sd(obsdifmet,na.rm=T)

	# calculate 2-D fit
	if(fitexpmeans){
                fit2d=lm(expmeanflux ~ expmeanconc)
        } else {
                fit2d=lm(yvar ~ xvar2d)
        }

        slope = fit2d$coef[2]; intercept = fit2d$coef[1] ## slope = 1 and intercept = 0 bc xvar generated from mlrfit
	print(paste("2D slope=", round(slope,2), "2D intercept",  round(intercept,2)))

	X_2 = seq(min(c(xvar2d,xlm))-0.1*diff(range(c(xvar2d,xlm))),max(c(xvar2d,xlm))+0.1*diff(range(c(xvar2d,xlm))),length.out=100) 
	y_fitmodel = X_2*slope+intercept

	beta=mvrnorm(n=nsamp2d,mu=as.vector(fit2d$coef),Sigma=vcov(fit2d))
	#mvrnorm(n = 1, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
	#Arguments:
	#       n: the number of samples required.
	#      mu: a vector giving the means of the variables.
	#   Sigma: a positive-definite symmetric matrix specifying the covariance matrix of the variables.

	# calc 2-D fit uncertainty using correlated sampling of slope and intercept errors
    	all_plot = matrix(0,nrow=nsamp2d, ncol=length(X_2))
	fitatobs=NULL
	for(s in c(1:nsamp2d)){
	        all_plot[s,] = X_2*beta[s,2]+beta[s,1]
		fitatobs=c(fitatobs,obsval2d*beta[s,2]+beta[s,1])

	}
	fitunc2d=sd(fitatobs)
	print(paste('Fit uncertainty at obs (2D) =',round(fitunc2d,3)))

	## to make flux estimate, need to also incorporate uncertainty in obs and PI
	## construct ensemble by
	# 1) calc PI from 2D fit
	# 2) for the primary sensor:
       		# picking one of 1,000 sets of coefficients to calculate predicted flux from Cset, Ct, and Cnet, and
		# picking one of 1,000 samples from prediction interval to offset (assuming PI about mean fit doesn't change significantly for different fits - would be more work to recalc PI every time)
	# 3) do the same for four simulated sensors constructed as primary+(sensor-NOAA)
	# 4) use primary for central estimate and SD of all (primary and simulated) for uncertainty estimate

	# 1) 
	y_err=sd(summary(fit2d)$res) # y component of the residuals

	# 2)
	allpredflux=NULL
	allpredfluxnopi=NULL
	beta=mvrnorm(n=nsampmlr,mu=as.vector(mlrfit$coef),Sigma=vcov(mlrfit))
	# first do for CO2_OP
	for(i in c(1:nsampmlr)){ # 2) nsampmlr random sets of fit coef
		for(j in c(1:nsampmlr)){ # 3) nsampmlr random samples of PI
			PI = randn()*y_err # get a new PI for every sensor and fit coef ## randn() has mean of 0 and SD of 1, and y_err is SD of residuals, so PI has a mean of 0 and an SD of y_err
			## so sampling from randn()*y_err should give the same spread as randomly picking one model's residual
			if(!is.null(xvar3)){
				allpredflux=c(allpredflux,beta[i,1] + beta[i,2] * obsval1[variables==varsel] + beta[i,3] * obsval2[variables==varsel] + beta[i,4] * obsval3[variables==varsel] + PI)
				allpredfluxnopi=c(allpredfluxnopi,beta[i,1] + beta[i,2] * obsval1[variables==varsel] + beta[i,3] * obsval2[variables==varsel] + beta[i,4] * obsval3[variables==varsel])
			} else if(!is.null(xvar2)){
				allpredflux=c(allpredflux,beta[i,1] + beta[i,2] * obsval1[variables==varsel] + beta[i,3] * obsval2[variables==varsel] + PI)
				allpredfluxnopi=c(allpredfluxnopi,beta[i,1] + beta[i,2] * obsval1[variables==varsel] + beta[i,3] * obsval2[variables==varsel])
			} else {
				allpredflux=c(allpredflux,beta[i,1] + beta[i,2] * obsval1[variables==varsel] + PI)
				allpredfluxnopi=c(allpredfluxnopi,beta[i,1] + beta[i,2] * obsval1[variables==varsel] + PI)
			}
		}
	}
	fluxnoaa=mean(allpredflux)
	fluxnoaasd=sd(allpredflux,na.rm=T)
	print(paste('*** CO2_OP only = ',round(fluxnoaa,2),' PgC/yr (+/- ',round(fluxnoaasd,2),')',sep=''))
	# 2a) for CO2_OP calc SD of just fit uncertainty at obs value - compare to SD of 2-fit uncertainty at obs line
	fituncmlr=sd(allpredfluxnopi)
	print(paste('Fit uncertainty at obs (MLR) =',round(fituncmlr,3)))

	# 3) next do for 4 simulated sensors (sensor-NOAA plus CO2_OP)
	for(difvar in sdsel){ # 1) pick 1 of 4 sets of obsvals (CO2_OP + 4 sensor-NOAA diffs)
		obsval1sim=obsval1[variables==varsel] + obsval1[variables==difvar]
		if(!is.null(xvar2)) obsval2sim=obsval2[variables==varsel] + obsval2[variables==difvar]
		if(!is.null(xvar3)) obsval3sim=obsval3[variables==varsel] + obsval3[variables==difvar]
		beta=mvrnorm(n=nsampmlr,mu=as.vector(mlrfit$coef),Sigma=vcov(mlrfit)) # get a new batch of nsampmlr for every sensor
		for(i in c(1:nsampmlr)){ # 2) nsampmlr random sets of fit coef 
			for(j in c(1:nsampmlr)){ # 3)  sampling PI
				PI = randn()*y_err # get a new PI for every sensor and fit coef
				if(!is.null(xvar3)){
					allpredflux=c(allpredflux,beta[i,1] + beta[i,2] * obsval1sim + beta[i,3] * obsval2sim + beta[i,4] * obsval3sim + PI)
				} else if(!is.null(xvar2)){
					allpredflux=c(allpredflux,beta[i,1] + beta[i,2] * obsval1sim + beta[i,3] * obsval2sim + PI)
				} else {
					allpredflux=c(allpredflux,beta[i,1] + beta[i,2] * obsval1sim + PI)
				}
			}
		}
	}

	# 4)
	fluxall5=mean(allpredflux) 
	fluxall5sd=sd(allpredflux,na.rm=T) 
	
	flux=fluxnoaa # using CO2_OP-only flux
	fluxsd=fluxall5sd # using sd of all

	print(paste('*** MLR Monte Carlo Flux = ',round(fluxall5,2),' PgC/yr (+/- ',round(fluxall5sd,2),')',sep=''))
	print(paste('Original Full SD =',round(sd(yvar),2)))
	print(paste('*** MLR Uncertainty Reduction =',round((1-fluxall5sd/sd(yvar))*100,1)))
	print(paste('Y_err =',round(y_err,2)))
	print(paste('SD of models within Obs 1SD',round(sd(yvar[xvar2d>(obsval2d-obssd)&xvar2d<(obsval2d+obssd)]),2)))

	plot_line_plus = rep(0,length(X_2))
	plot_line_minus = rep(0,length(X_2))
	for(i in c(1:length(X_2))){
		plot_line_plus[i] = mean(all_plot[,i]) + (var(all_plot[,i]) + y_err^2)^0.5
		plot_line_minus[i] = mean(all_plot[,i]) - (var(all_plot[,i]) + y_err^2)^0.5
	}

	library('RColorBrewer')
	cols=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2')) # 17 diff colors
	#cols[6]='gold' # yellow too hard to see

	experiments=c('IS','LNLG','LNLGIS','LNLGOGIS','OG') # ordered to minimize lines doubling back
	expcols=cols[c(1,3,5,4,2)]
	expcols=expcols[is.element(experiments,expnames)]
	experiments=experiments[is.element(experiments,expnames)]
	models=unique(modnames)


	pngdim=1200
	pngpts=24

	for(i in c(1:length(plotsel))){ # now only allowing single panel plots but allowing multiple in one call

		filename=filenames[i]

		png(filename,height=pngdim,width=pngdim,pointsize=pngpts) # *length(plotsel),pointsize=pngpts)
		par(mfrow=c(1,1)) # length(plotsel)))
		par(mgp=c(2.75,1,0))
		par(mar=c(5,4.5,2,0)+0.1)
		par(oma=c(0,0,3,0))

		if(is.null(xlm)) xlm=range(c(xvar2d,obsval2d-obssd,obsval2d+obssd))
		if(is.null(ylm)) ylm=range(yvar)

		if(plotsel[i]==1){

			plot(xvar2d,yvar,main='By Experiment',type='n',ylab=ylb,xlab=xlb,cex.main=2,cex.axis=1.5,cex.lab=1.5,xlim=xlm,ylim=ylm)
			y1=par('usr')[3]; y2=par('usr')[4]
			polygon(c(obsval2d-obssd,obsval2d+obssd,obsval2d+obssd,obsval2d-obssd,obsval2d-obssd),c(y1,y1,y2,y2,y1),col='gray80',border=NA)
			abline(v=obsval2d,col='black',lwd=3)
			for(mod in models){
				points(xvar2d[modnames==mod],yvar[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
			}
			abline(intercept,slope,col='Black',lwd=6)
			points(expmeanconc,expmeanflux,cex=2,col=expcols[1:length(experiments)],pch=16)
			segments(expmeanconc-expsdconc,expmeanflux,expmeanconc+expsdconc,expmeanflux,lwd=4,col=expcols[1:length(experiments)])
			segments(expmeanconc,expmeanflux-expsdflux,expmeanconc,expmeanflux+expsdflux,lwd=4,col=expcols[1:length(experiments)])
		#	mtext(paste('m =',round(slope,2)),1,-3,col='Black')
		#	mtext(paste('r^2 = ',round(r2,2),' (',round(cor(expmeanconc,expmeanflux)^2,2),' for exp means)',sep=''),1,-4.5)
		# 	mtext(paste('uncred =',ur),1,-1.5)
			if(inclmodleg) legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
			if(inclexpleg) legend('bottomright',c(experiments,'fit'),col=c(expcols[c(1:length(experiments))],'Black'),pt.cex=2.0,pt.lwd=3,pch=c(rep(16,length(experiments)),NA),pt.bg='white',bty='n',lwd=c(rep(F,length(experiments)),6),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
			if(!is.null(r2loc)) legend(r2loc,legend=substitute(paste(r^2,' = ',v),list(v=round(cor(xvar2d,yvar)^2,2))),bty='n',cex=1.5)
			print(paste('r^2 = ',round(cor(xvar2d,yvar)^2,3)))
			print(paste('exp mean r^2 =',round(cor(expmeanconc,expmeanflux)^2,3)))

		}
		if(plotsel[i]==2){

			plot(xvar2d,yvar,main='By Model',type='n',ylab=ylb,xlab=xlb,cex.main=2,cex.axis=1.5,cex.lab=1.5,xlim=xlm,ylim=ylm)
			y1=par('usr')[3]; y2=par('usr')[4]
			for(mod in models){
				points(xvar2d[modnames==mod],yvar[modnames==mod],col=cols[match(mod,models)],cex=2,lwd=3,pch=match(expnames[modnames==mod],experiments),bg='white',type='b')
			}
			abline(intercept,slope,col='Black',lwd=6)
			modmeanconc=aggregate(xvar2d,by=list(mod=modnames),mean)$x # alphabetizes model names
			modmeanflux=aggregate(yvar,by=list(mod=modnames),mean)$x # alphabetizes model names
			points(modmeanconc,modmeanflux,cex=2,col=cols[order(models)],pch=16)
			if(inclmodleg) legend('topleft',models,col=cols[c(1:length(models))],pt.cex=2.0,pt.lwd=3,pch=16,pt.bg='white',bty='n',lty=1,ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
			if(inclexpleg) legend('bottomright',c(experiments,'fit'),pch=c(c(1:length(experiments)),NA),col=c(rep('black',length(experiments)),'Black'),pt.cex=2.0,pt.lwd=3,pt.bg='white',bty='n',lwd=c(rep(F,length(experiments)),6),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
			if(!is.null(r2loc)) legend(r2loc,legend=substitute(paste(r^2,' = ',v),list(v=round(cor(xvar2d,yvar)^2,2))),bty='n',cex=1.5)
			box()

		}
		if(plotsel[i]==3){

			plot(xvar2d,yvar,main=ecmain,type='n',ylab=ylb,xlab=xlb,cex.main=2,cex.axis=1.5,cex.lab=1.5,xlim=xlm,ylim=ylm)
			y1=par('usr')[3]; y2=par('usr')[4]
			polygon(c(obsval2d-obssd,obsval2d+obssd,obsval2d+obssd,obsval2d-obssd,obsval2d-obssd),c(y1,y1,y2,y2,y1),col='gray80',border=NA)
			abline(v=obsval2d,col='black',lwd=3)
			polygon(c(X_2,rev(X_2)),c(plot_line_plus,rev(plot_line_minus)),col='gold',border=NA,density=50)
			for(mod in models){
				points(xvar2d[modnames==mod],yvar[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
			}
			abline(intercept,slope,col='Black',lwd=6)
			points(expmeanconc,expmeanflux,cex=2,col=expcols[1:length(experiments)],pch=16)
			segments(expmeanconc-expsdconc,expmeanflux,expmeanconc+expsdconc,expmeanflux,lwd=4,col=expcols[1:length(experiments)])
			segments(expmeanconc,expmeanflux-expsdflux,expmeanconc,expmeanflux+expsdflux,lwd=4,col=expcols[1:length(experiments)])
			if(inclmodleg) legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
			if(inclexpleg) legend('bottomright',c(experiments,'fit'),col=c(expcols[c(1:length(experiments))],'Black'),pt.cex=2.0,pt.lwd=3,pch=c(rep(16,length(experiments)),NA),pt.bg='white',bty='n',lwd=c(rep(F,length(experiments)),6),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0.1)
			if(!is.null(r2loc)) legend(r2loc,legend=substitute(paste(r^2,' = ',v),list(v=round(cor(xvar2d,yvar)^2,2))),bty='n',cex=1.5)
			box()

			## EPS
			width_cm <- 8.9
			width_in <- width_cm / 2.54
			sc <- 0.21 # Scaling factor for absolute line widths
			new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)
			file_name=paste(substr(filename,1,nchar(filename)-3),'eps',sep='')
			cairo_ps(filename = file_name,
			         width = width_in,
			         height = width_in,
			         pointsize = new_pointsize)

			plot(xvar2d,yvar,main=ecmain,type='n',ylab='',xlab=xlb,cex.main=2,cex.axis=1.5,cex.lab=1.5,xlim=xlm,ylim=ylm,axes=F)
			mtext(ylb,2,2.3,cex=1.5)
			axis(1,cex.axis=1.5,lwd=0.5)
			axis(2,cex.axis=1.5,lwd=0.5)

                        y1=par('usr')[3]; y2=par('usr')[4]
                        polygon(c(obsval2d-obssd,obsval2d+obssd,obsval2d+obssd,obsval2d-obssd,obsval2d-obssd),c(y1,y1,y2,y2,y1),col='gray80',border=NA)
                        abline(v=obsval2d,col='black',lwd=3*sc)
                        polygon(c(X_2,rev(X_2)),c(plot_line_plus,rev(plot_line_minus)),col=lighten('gold',0.25),border=NA) #,density=50)
                        for(mod in models){
                                points(xvar2d[modnames==mod],yvar[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3*sc,pch=match(mod,models),bg='white')
                        }
                        abline(intercept,slope,col='Black',lwd=6*sc)
                        points(expmeanconc,expmeanflux,cex=2,col=expcols[1:length(experiments)],pch=16)
                        segments(expmeanconc-expsdconc,expmeanflux,expmeanconc+expsdconc,expmeanflux,lwd=4*sc,col=expcols[1:length(experiments)])
                        segments(expmeanconc,expmeanflux-expsdflux,expmeanconc,expmeanflux+expsdflux,lwd=4*sc,col=expcols[1:length(experiments)])
                        if(inclmodleg) legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3*sc,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
                        if(inclexpleg) legend('bottomright',c(experiments,'fit'),col=c(expcols[c(1:length(experiments))],'Black'),pt.cex=2.0,pt.lwd=3*sc,pch=c(rep(16,length(experiments)),NA),pt.bg='white',bty='n',lwd=c(rep(NA,length(experiments)),6*sc),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0.1)
                        if(!is.null(r2loc)) legend(r2loc,legend=substitute(paste(r^2,' = ',v),list(v=round(cor(xvar2d,yvar)^2,2))),bty='n',cex=1.5)
                        box(lwd=0.5)
			dev.off()

		}
		if(plotsel[i]==4){

			plot(xvar2d,yvar,main='By Experiment',type='n',ylab=ylb,xlab=xlb,cex.main=2,cex.axis=1.5,cex.lab=1.5,xlim=xlm,ylim=ylm)
			y1=par('usr')[3]; y2=par('usr')[4]
			for(mod in models){
				points(xvar2d[modnames==mod],yvar[modnames==mod],col=expcols[match(expnames[modnames==mod],experiments)],cex=2,lwd=3,pch=match(mod,models),bg='white')
			}
			abline(intercept,slope,col='Black',lwd=6)
			points(expmeanconc,expmeanflux,cex=2,col=expcols[1:length(experiments)],pch=16)
			segments(expmeanconc-expsdconc,expmeanflux,expmeanconc+expsdconc,expmeanflux,lwd=4,col=expcols[1:length(experiments)])
			segments(expmeanconc,expmeanflux-expsdflux,expmeanconc,expmeanflux+expsdflux,lwd=4,col=expcols[1:length(experiments)])
			if(inclmodleg) legend('topleft',models,col='black',pt.cex=2.0,pt.lwd=3,pch=c(1:length(models)),pt.bg='white',bty='n',ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
			if(inclexpleg) legend('bottomright',c(experiments,'fit'),col=c(expcols[c(1:length(experiments))],'Black'),pt.cex=1.5,pt.lwd=3,pch=c(rep(16,length(experiments)),NA),pt.bg='white',bty='n',lwd=c(rep(F,length(experiments)),6),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0.1)
			if(!is.null(r2loc)) legend(r2loc,legend=substitute(paste(r^2,' = ',v),list(v=round(cor(xvar2d,yvar)^2,2))),bty='n',cex=1.5)

		}
		if(plotsel[i]==5){

			plot(xvar2d,yvar,main=ecmain,type='n',ylab=ylb,xlab=xlb,cex.main=2,cex.axis=1.5,cex.lab=1.5,xlim=xlm,ylim=ylm)
			y1=par('usr')[3]; y2=par('usr')[4]
			polygon(c(obsval2d-obssd,obsval2d+obssd,obsval2d+obssd,obsval2d-obssd,obsval2d-obssd),c(y1,y1,y2,y2,y1),col='gray80',border=NA)
			abline(v=obsval2d,col='black',lwd=3)
			polygon(c(X_2,rev(X_2)),c(plot_line_plus,rev(plot_line_minus)),col='gold',border=NA,density=50)
			# from ATOM/SCI/CASEAS/MODELS/OCO2V10MIP/oco_flux_vs_flux.r
			#modord=order(modsdland,decreasing=T)
			#for(mod in mods[modord]){
			for(mod in models){
				expord=order(xvar2d[modnames==mod])
				points(xvar2d[modnames==mod][expord],yvar[modnames==mod][expord],col=cols[match(mod,models)],cex=2,lwd=3,pch=match(expnames[modnames==mod][expord],experiments),bg='white',type='b')
			}
			abline(intercept,slope,col='Black',lwd=6)
			modmeanconc=aggregate(xvar2d,by=list(mod=modnames),mean)$x # alphabetizes model names
			modmeanflux=aggregate(yvar,by=list(mod=modnames),mean)$x # alphabetizes model names
	#                points(modmeanconc,modmeanflux,cex=2,col=cols[order(models)],pch=16)
			if(inclmodleg) legend('topleft',models,col=cols[c(1:length(models))],pt.cex=2.0,pt.lwd=3,pch=16,pt.bg='white',bty='n',lty=1,ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
			if(inclexpleg) legend('bottomright',c(experiments,'fit'),pch=c(c(1:length(experiments)),NA),col=c(rep('black',length(experiments)),'Black'),pt.cex=1.5,pt.lwd=3,pt.bg='white',bty='n',lwd=c(rep(F,length(experiments)),6),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
			if(!is.null(r2loc)) legend(r2loc,legend=substitute(paste(r^2,' = ',v),list(v=round(cor(xvar2d,yvar)^2,2))),bty='n',cex=1.5)
			box()

			## EPS
                        width_cm <- 8.9
                        width_in <- width_cm / 2.54
                        sc <- 0.21 # Scaling factor for absolute line widths
                        new_pointsize <- 24 * sc # ~5.04 (This handles the base text size)
                        file_name=paste(substr(filename,1,nchar(filename)-3),'eps',sep='')
                        cairo_ps(filename = file_name,
                                 width = width_in,
                                 height = width_in,
                                 pointsize = new_pointsize)
			plot(xvar2d,yvar,main=ecmain,type='n',ylab='',xlab=xlb,cex.main=2,cex.axis=1.5,cex.lab=1.5,xlim=xlm,ylim=ylm,axes=F)
			mtext(ylb,2,2.3,cex=1.5)
                        axis(1,cex.axis=1.5,lwd=0.5)
                        axis(2,cex.axis=1.5,lwd=0.5)
                        y1=par('usr')[3]; y2=par('usr')[4]
                        polygon(c(obsval2d-obssd,obsval2d+obssd,obsval2d+obssd,obsval2d-obssd,obsval2d-obssd),c(y1,y1,y2,y2,y1),col='gray80',border=NA)
                        abline(v=obsval2d,col='black',lwd=3*sc)
                        polygon(c(X_2,rev(X_2)),c(plot_line_plus,rev(plot_line_minus)),col=lighten('gold',0.25),border=NA) # ,density=50)
                        for(mod in models){
                                expord=order(xvar2d[modnames==mod])
                                points(xvar2d[modnames==mod][expord],yvar[modnames==mod][expord],col=cols[match(mod,models)],cex=2,lwd=3*sc,pch=match(expnames[modnames==mod][expord],experiments),bg='white',type='b')
                        }
                        abline(intercept,slope,col='Black',lwd=6*sc)
                        modmeanconc=aggregate(xvar2d,by=list(mod=modnames),mean)$x # alphabetizes model names
                        modmeanflux=aggregate(yvar,by=list(mod=modnames),mean)$x # alphabetizes model names
                        if(inclmodleg) legend('topleft',models,col=cols[c(1:length(models))],pt.cex=2.0,pt.lwd=3*sc,pch=16,pt.bg='white',bty='n',lty=1,lwd=3*sc,ncol=1,bg='white',box.col='white',y.intersp=1.25,title='Model:',title.col='black',title.adj=0,inset=0.02)
                        if(inclexpleg) legend('bottomright',c(experiments,'fit'),pch=c(c(1:length(experiments)),NA),col=c(rep('black',length(experiments)),'Black'),pt.cex=1.5,pt.lwd=3*sc,pt.bg='white',bty='n',lwd=c(rep(NA,length(experiments)),6*sc),bg='white',box.col='white',title='Experiment:',title.col='black',title.adj=0)
                        if(!is.null(r2loc)) legend(r2loc,legend=substitute(paste(r^2,' = ',v),list(v=round(cor(xvar2d,yvar)^2,2))),bty='n',cex=1.5)
                        box(lwd=0.5)
			dev.off()

		}

		dev.off()
	
	} # loop through plotsel

	return(c(flux,fluxsd))

} # end of function
