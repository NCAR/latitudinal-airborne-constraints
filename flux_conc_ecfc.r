## program to correlate model posterior CO2 and fluxes for large press x lat regions for 50 models against 50 mean fluxes and plot and output ECFC results

## must run flux_conc_cor.r first


library(abind)
library(pracma)
source('ecplot_mlr.r')

## Options:
fitmodmeans=F; fitexpmeans=F # options instead of individual inversions, can not both be true

exclog=F # whether to exclude OG and LNLGOGIS
oneexponly=NULL # 'IS'  
rsqco=0

latlimfluxes=F # whether to trim fluxes to latvals (from BINFIT concentration output)

# the 2D fit is sampled just nsample2d times, so not time consuming:
nsample2d=1000

# the MLR fits get sampled 5*nsamplemlr^2 times, so time consuming:
# nsamplemlr=1000 stalls on first set (SET) (100% CPU at 12 hours)

# takes 384 sec:
nsamplemlr=200 ## default for review response
#nsamplemlr=50 # for testing

# for the diff metrics (e.g. SET-T conc vs SET-T flux), using fewer iterations:
#nsamplemlrfast=100
nsamplemlrfast=50


models=c('AMES','BAKER','CAMS','CMS-FLUX','COLA','CT','OU','TM5-4DVAR','UT','WOMBAT')

modeldir='V10MIP'
obsdir='OBS/BINFIT'

## set plot size and axis ranges
pngw=650;pngh=650;psize=12
xlim1=c(-70,80);ylim1=c(1000,300)
xtick=20

## read in results of stacking 2D (P vs lat) annual mean fields for each model in flux_conc_cor.r

temp=read.table(paste(modeldir,'/Annual_mean_2D_fields_as_columns.txt',sep=''),header=T)
allannmat=temp[,3:ncol(temp)] # 2D fields expanded into single columns
modnames=gsub('\\.','-',sapply(strsplit(colnames(temp[3:ncol(temp)]),'_'),"[[", 1))
expnames=sapply(strsplit(colnames(temp[3:ncol(temp)]),'_'),"[[", 2)
latvals=scan(paste(modeldir,'/',models[1],'/IS/Annual_Mean_CO2_OP_10x100.txt',sep=''),nlines=1) # same for all experiments
prsvals=scan(paste(modeldir,'/',models[1],'/IS/Annual_Mean_CO2_OP_10x100.txt',sep=''),skip=1,nlines=1)
alllat=temp[,1]
allprs=temp[,2]
annbg=scan(paste(modeldir,'/Annual_mean_2D_field_averages.txt',sep=''))

temp=read.table(paste(obsdir,'/Annual_mean_2D_fields_as_columns.txt',sep=''),header=T)
obsannmat=temp[,3:ncol(temp)] # 2D fields expanded into single columns
variables=colnames(temp[3:ncol(temp)])


## set plot colors
mattcol=read.table('/h/eol/stephens/NCAR/ORCAS/PROC/MERGE/mattcol127.txt',header=F)
colnames(mattcol)=c('red','green','blue')
colormc=rgb(mattcol,maxColorValue=255)

## read in fluxes output by V10MIP/flux_v_lat.r
fluxes=read.table(paste(modeldir,'/oco2_v10mip_total_fluxes_bylat_201606-201805.txt',sep=''),header=T,stringsAsFactors=F)

if(exclog){
	fluxes=fluxes[!grepl('OG',fluxes$exp),]
	allannmat=allannmat[,expnames!='OG']
	modnames=modnames[expnames!='OG']
	expnames=expnames[expnames!='OG']
}
if(!is.null(oneexponly)){
	fluxes=fluxes[grepl(oneexponly,fluxes$exp),]
	allannmat=allannmat[,expnames==oneexponly]
	modnames=modnames[expnames==oneexponly]
	expnames=expnames[expnames==oneexponly]
}

fluxes=fluxes[which(is.element(fluxes[,2],modnames)),] # gets rid of NIES, CSU, and WEIR
fluxes=fluxes[which(is.element(fluxes[,1],expnames)),] # gets rid of extra experiments if !is.null(oneexponly)

if(any(modnames!=fluxes[,2])){ stop('MISMATCH BETWEEN CONC AND FLUX ROWS') }
if(any(expnames!=fluxes[,1])){ stop('MISMATCH BETWEEN CONC AND FLUX ROWS') }

fluxlats=substr(names(fluxes),nchar(names(fluxes))-2,nchar(names(fluxes)))
fluxlats=fluxlats[3:length(fluxlats)]
fluxlats=gsub('\\.','-',fluxlats)
fluxlats=gsub('X','',fluxlats)
fluxlats=as.numeric(fluxlats)

fluxnums=t(fluxes[,3:ncol(fluxes)])

if(latlimfluxes){

	## binfit was done from 70S to 80N but fluxes go from 90S to 90N
	fluxnums=fluxnums[fluxlats>=min(latvals)&fluxlats<=max(latvals),]
	fluxlats=fluxlats[fluxlats>=min(latvals)&fluxlats<=max(latvals)]

}

## calc fluxes for 6 regions = SET, T, NET, SET-T, T-(SET+NET), NET-T
setflux=apply(fluxnums[fluxlats<=(-20),],2,sum)
tropflux=apply(fluxnums[fluxlats>=(-20)&fluxlats<=20,],2,sum)
netflux=apply(fluxnums[fluxlats>=20,],2,sum)
setmtropflux=setflux-tropflux
tropmsetnetflux=tropflux-(setflux+netflux)
netmtropflux=netflux-tropflux
regfluxes=cbind(setflux,tropflux,netflux,setmtropflux,tropmsetnetflux,netmtropflux)
labels=c('SET','Tropical','NET','SET-Trop','Trop-(SET+NET)','NET-Trop')
labels2=format(labels,width=max(nchar(labels)))

## calc annmn fits for each bin against fluxes from 6 regions
rsqds=matrix(NA,ncol=6,nrow=nrow(allannmat))
coefs=matrix(NA,ncol=6,nrow=nrow(allannmat))
for(i in c(1:nrow(allannmat))){ # loop on conc bins
	for(j in c(1:6)){
		fit=lm(unlist(allannmat[i,]) ~ regfluxes[,j])
		rsqds[i,j]=summary(fit)$r.squared
		coefs[i,j]=fit$coef[2]
	}
}


## Calculate regional mean concentrations

scalefactor=rep(cos(latvals/180*pi),length(prsvals)) # converting to PgC equiv

setsel<-alllat<=-30
tropsel<-alllat>=-10&alllat<=0
netsel<-alllat>=20
etsel<-setsel|netsel
latdivs='-30, -10 to 0, 20'

# in order to apply rsqco, need to calc for each flux region separately (as above)
annset=NULL
anntrop=NULL
annnet=NULL
annet=NULL
# allannmat is 105 bins x 50 models
for(i in c(1:6)){
        rsqdsel<-rsqds[,i]>=rsqco
        if(any(rsqdsel)){
                tmpset=apply(rsqdsel*setsel*scalefactor*allannmat,2,sum)/sum(rsqdsel*setsel*scalefactor)
                tmptrop=apply(rsqdsel*tropsel*scalefactor*allannmat,2,sum)/sum(rsqdsel*tropsel*scalefactor)
                tmpnet=apply(rsqdsel*netsel*scalefactor*allannmat,2,sum)/sum(rsqdsel*netsel*scalefactor)
                tmpet=apply(rsqdsel*etsel*scalefactor*allannmat,2,sum)/sum(rsqdsel*etsel*scalefactor)
                annset=cbind(annset,tmpset)
                anntrop=cbind(anntrop,tmptrop) ## e.g. anntrop will be different for set-trop and for net-trop bc rsqco cuts differently
                annnet=cbind(annnet,tmpnet)
                annet=cbind(annet,tmpet)
        }
}

## do for obs too

obsannset=NULL
obsanntrop=NULL
obsannnet=NULL
obsannet=NULL
obsannglb=NULL
# obsannmat is 105 bins x 5 variables
for(i in c(1:6)){ ## 6 would only be different if rsqco!=0
	rsqdsel<-rsqds[,i]>=rsqco
        if(any(rsqdsel)){
                tmpset=apply(rsqdsel*setsel*scalefactor*obsannmat,2,sum,na.rm=T)/sum(rsqdsel*setsel*scalefactor)
                tmptrop=apply(rsqdsel*tropsel*scalefactor*obsannmat,2,sum,na.rm=T)/sum(rsqdsel*tropsel*scalefactor)
                tmpnet=apply(rsqdsel*netsel*scalefactor*obsannmat,2,sum,na.rm=T)/sum(rsqdsel*netsel*scalefactor)
                tmpet=apply(rsqdsel*etsel*scalefactor*obsannmat,2,sum,na.rm=T)/sum(rsqdsel*etsel*scalefactor)
                tmpglb=apply(rsqdsel*scalefactor*obsannmat,2,sum,na.rm=T)/sum(rsqdsel*scalefactor)
                obsannset=cbind(obsannset,tmpset)
                obsanntrop=cbind(obsanntrop,tmptrop)
                obsannnet=cbind(obsannnet,tmpnet)
                obsannet=cbind(obsannet,tmpet)
                obsannglb=cbind(obsannglb,tmpglb)
        }
}

## for obs and obs error:
sdsel=c('CO2_AO2_m_CO2_NOAA','CO2_QCLS_m_CO2_NOAA','CO2_MED_m_CO2_NOAA','CO2_PFP_m_CO2_NOAA')
obsdifremppm=NULL # c(sdsel,'CO2_NOAA_m_CO2_NOAA')
obsdifrempgc=NULL # c(sdsel,'CO2_NOAA_m_CO2_NOAA')

print(paste('(assuming rsqco=0) Observed SET, T, NET, ET =',paste(round(obsannset[variables=='CO2_OP',1],3),round(obsanntrop[variables=='CO2_OP',1],3),round(obsannnet[variables=='CO2_OP',1],3),round(obsannet[variables=='CO2_OP',1],3),collapse=', ')))
obsconc=c(round(obsannset[variables=='CO2_OP',1],3),round(obsanntrop[variables=='CO2_OP',1],3),round(obsannnet[variables=='CO2_OP',1],3),round(obsannet[variables=='CO2_OP',1],3),round(obsannglb[variables=='CO2_OP',1],3))
## if rsqco!=0 then obs value will change depending on which flux it is being correlated to 


print('SET')
print('2-Var MLR')
i=1
fit=lm(regfluxes[,i] ~ annset[,i] + anntrop[,i])
print(paste('Adjusted R-squared =',round(summary(fit)$adj.r.squared,3)))
# calc experiment means using only models with all experiments present
expcount=aggregate(!is.na(regfluxes[,i]),by=list(mod=modnames),sum)
sel<-expcount$x[match(modnames,expcount$mod)]==length(unique(expnames))
expmeanannset=aggregate(annset[,i][sel],by=list(exp=expnames[sel]),mean)$x
expmeananntrop=aggregate(anntrop[,i][sel],by=list(exp=expnames[sel]),mean)$x
expmeanflux=aggregate(regfluxes[,i][sel],by=list(exp=expnames[sel]),mean)$x
print(paste('Exp Mean Adjusted R-squared =',round(summary(lm(expmeanflux ~ expmeanannset + expmeananntrop))$adj.r.squared,3)))
print(summary(fit)$coef)
print(paste(labels2[i],round(summary(fit)$adj.r.squared,3),round((1-sd(fit$resid)/sd(regfluxes[,i]))*100,2),round((1-sd(fit$resid)/sd(regfluxes[!grepl('OG',expnames),i]))*100,2)))
obsdifmet=NULL
for(j in c(1:length(sdsel))){
        obsdifmet=c(obsdifmet,(obsannset[variables==sdsel[j],i]*fit$coef[2]+obsanntrop[variables==sdsel[j],i]*fit$coef[3])) ## not adding in fit$coef[1] because only interested in difference (would need to use fit$coef[1] instead of 0 for NOAA-NOAA if we did)
}
obsdifmet=c(obsdifmet,0) # add 0 for NOAA-NOAA
obsdifsd=sd(obsdifmet,na.rm=T)
names(obsdifmet)=c(sdsel,'CO2_NOAA_m_CO2_NOAA')
print(obsdifmet)
print(obsdifsd)
obsdifrempgc=cbind(obsdifrempgc,round(obsdifmet,3))
table2remcoef=c(fit$coef[2],fit$coef[3],NA,fit$coef[1])
table2remse=c(summary(fit)$coef[2,2],summary(fit)$coef[3,2],NA,summary(fit)$coef[1,2])

fluxout=ecplotmlr(xvar1=annset[,i],xvar2=anntrop[,i],xvar3=NULL,yvar=regfluxes[,i],obsval1=obsannset[,i],obsval2=obsanntrop[,i],obsval3=NULL,nsamp2d=nsample2d,nsampmlr=nsamplemlr,modnames=modnames,expnames=expnames,filenames=c(paste(modeldir,'/SET_Flux_vs_SET-T_Conc_latdiv_MLR_fig4a.png',sep=''),paste(modeldir,'/SET_Flux_vs_SET-T_Conc_latdiv_MLR_figS4a.png',sep='')),plotsel=c(3,5),fitexpmeans=F,ecmain=expression('a. Southern Extratropics'),xlb=expression(paste('Concentration-based Flux Prediction (PgC ',yr^-1,')')),xlm=range(regfluxes[,i]),ylm=range(regfluxes[,i]),inclmodleg=F,inclexpleg=F,ylb=expression(paste('v10 MIP Total Flux (PgC ',yr^-1,')')) )
remsetlatdiv=c(mean(regfluxes[,i]),sd(regfluxes[,i]),mean(regfluxes[!grepl('OG',expnames),i]),sd(regfluxes[!grepl('OG',expnames),i]),fluxout[1],fluxout[2],fit$coef[1],fit$coef[2],fit$coef[3],NA)
system(paste('cp ',modeldir,'/SET_Flux_vs_SET-T_Conc_latdiv_MLR_fig4a.eps ',modeldir,'/fig4a.eps',sep=''))
system(paste('cp ',modeldir,'/SET_Flux_vs_SET-T_Conc_latdiv_MLR_figS4a.eps ',modeldir,'/figS4a.eps',sep=''))

# calc simple diff
obsdifmet11=NULL
for(j in c(1:length(sdsel))){
        obsdifmet11=c(obsdifmet11,obsannset[variables==sdsel[j],i]-obsanntrop[variables==sdsel[j],i])
}
obsdifmet11=c(obsdifmet11,0) # add 0 for NOAA-NOAA
names(obsdifmet11)=c(sdsel,'CO2_NOAA_m_CO2_NOAA')
obsdifremppm=cbind(obsdifremppm,round(obsdifmet11,3))


print('Tropical')
print('3-Var MLR')
i=2
fit=lm(regfluxes[,i] ~ anntrop[,i] + annset[,i] + annnet[,i])
print(paste('Adjusted R-squared =',round(summary(fit)$adj.r.squared,3)))
# calc experiment means using only models with all experiments present
expcount=aggregate(!is.na(regfluxes[,i]),by=list(mod=modnames),sum)
sel<-expcount$x[match(modnames,expcount$mod)]==length(unique(expnames))
expmeanannnet=aggregate(annnet[,i][sel],by=list(exp=expnames[sel]),mean)$x
expmeanannset=aggregate(annset[,i][sel],by=list(exp=expnames[sel]),mean)$x
expmeananntrop=aggregate(anntrop[,i][sel],by=list(exp=expnames[sel]),mean)$x
expmeanflux=aggregate(regfluxes[,i][sel],by=list(exp=expnames[sel]),mean)$x
print(paste('Exp Mean Adjusted R-squared =',round(summary(lm(expmeanflux ~ expmeananntrop + expmeanannnet + expmeanannnet))$adj.r.squared,3)))
print(summary(fit)$coef)
print(paste(labels2[i],round(summary(fit)$adj.r.squared,3),round((1-sd(fit$resid)/sd(regfluxes[,i]))*100,2),round((1-sd(fit$resid)/sd(regfluxes[!grepl('OG',expnames),i]))*100,2)))
# 3 var:
obsdifmet=NULL
for(j in c(1:length(sdsel))){
        obsdifmet=c(obsdifmet,(obsanntrop[variables==sdsel[j],i]*fit$coef[2]+obsannset[variables==sdsel[j],i]*fit$coef[3]+obsannnet[variables==sdsel[j],i]*fit$coef[4])) # +fit$coef[1])) # )/mean(abs(fit$coef[2:3])))
}
obsdifmet=c(obsdifmet,0) # add 0 for NOAA-NOAA
obsdifsd=sd(obsdifmet,na.rm=T)
names(obsdifmet)=c(sdsel,'CO2_NOAA_m_CO2_NOAA')
print(obsdifmet)
print(obsdifsd)
obsdifrempgc=cbind(obsdifrempgc,round(obsdifmet,3))
table2remcoef=rbind(table2remcoef,c(fit$coef[3],fit$coef[2],fit$coef[4],fit$coef[1]))
table2remse=rbind(table2remse,c(summary(fit)$coef[3,2],summary(fit)$coef[2,2],summary(fit)$coef[4,2],summary(fit)$coef[1,2]))

fluxout=ecplotmlr(xvar1=anntrop[,i],xvar2=annset[,i],xvar3=annnet[,i],yvar=regfluxes[,i],obsval1=obsanntrop[,i],obsval2=obsannset[,i],obsval3=obsannnet[,i],nsamp2d=nsample2d,nsampmlr=nsamplemlr,modnames=modnames,expnames=expnames,filenames=c(paste(modeldir,'/T_Flux_vs_T-ET_Conc_latdiv_MLR_fig4b.png',sep=''),paste(modeldir,'/T_Flux_vs_T-ET_Conc_latdiv_MLR_figS4b.png',sep='')),plotsel=c(3,5),fitexpmeans=F,inclmodleg=T,inclexpleg=T,ecmain=expression('b. Tropics'),xlb=expression(paste('Concentration-based Flux Prediction (PgC ',yr^-1,')')),ylb='',xlm=range(regfluxes[,i])+c(0,0.25),ylm=range(regfluxes[,i])+c(0,0.25) ) # ,ylb=expression(paste('v10 MIP Total Flux (PgC ',yr^-1,')')) )
remtlatdiv=c(mean(regfluxes[,i]),sd(regfluxes[,i]),mean(regfluxes[!grepl('OG',expnames),i]),sd(regfluxes[!grepl('OG',expnames),i]),fluxout[1],fluxout[2],fit$coef[1],fit$coef[2],fit$coef[3],fit$coef[4])
system(paste('cp ',modeldir,'/T_Flux_vs_T-ET_Conc_latdiv_MLR_fig4b.eps ',modeldir,'/fig4b.eps',sep=''))
system(paste('cp ',modeldir,'/T_Flux_vs_T-ET_Conc_latdiv_MLR_figS4b.eps ',modeldir,'/figS4b.eps',sep=''))

# calc simple diff
obsdifmet11=NULL
for(j in c(1:length(sdsel))){
        obsdifmet11=c(obsdifmet11,obsanntrop[variables==sdsel[j],i]-obsannet[variables==sdsel[j],i])
	## averaging set and net together into et
}
obsdifmet11=c(obsdifmet11,0) # add 0 for NOAA-NOAA
names(obsdifmet11)=c(sdsel,'CO2_NOAA_m_CO2_NOAA')
obsdifremppm=cbind(obsdifremppm,round(obsdifmet11,3))


print('NET')
print('2-Var MLR')
i=3
fit=lm(regfluxes[,i] ~ annnet[,i] + anntrop[,i])
print(paste('Adjusted R-squared =',round(summary(fit)$adj.r.squared,3)))
# calc experiment means using only models with all experiments present
expcount=aggregate(!is.na(regfluxes[,i]),by=list(mod=modnames),sum)
sel<-expcount$x[match(modnames,expcount$mod)]==length(unique(expnames))
expmeanannnet=aggregate(annnet[,i][sel],by=list(exp=expnames[sel]),mean)$x
expmeananntrop=aggregate(anntrop[,i][sel],by=list(exp=expnames[sel]),mean)$x
expmeanflux=aggregate(regfluxes[,i][sel],by=list(exp=expnames[sel]),mean)$x
print(paste('Exp Mean Adjusted R-squared =',round(summary(lm(expmeanflux ~ expmeanannnet + expmeananntrop))$adj.r.squared,3)))
print(summary(fit)$coef)
print(paste(labels2[i],round(summary(fit)$adj.r.squared,3),round((1-sd(fit$resid)/sd(regfluxes[,i]))*100,2),round((1-sd(fit$resid)/sd(regfluxes[!grepl('OG',expnames),i]))*100,2)))
obsdifmet=NULL
for(j in c(1:length(sdsel))){
        obsdifmet=c(obsdifmet,(obsannnet[variables==sdsel[j],i]*fit$coef[2]+obsanntrop[variables==sdsel[j],i]*fit$coef[3])) # +fit$coef[1])) # )/mean(abs(fit$coef[2:3])))
}
obsdifmet=c(obsdifmet,0) # add 0 for NOAA-NOAA
obsdifsd=sd(obsdifmet,na.rm=T)
names(obsdifmet)=c(sdsel,'CO2_NOAA_m_CO2_NOAA')
print(obsdifmet)
print(obsdifsd)
obsdifrempgc=cbind(obsdifrempgc,round(obsdifmet,3))
table2remcoef=rbind(table2remcoef,c(NA,fit$coef[3],fit$coef[2],fit$coef[1]))
table2remse=rbind(table2remse,c(NA,summary(fit)$coef[3,2],summary(fit)$coef[2,2],summary(fit)$coef[1,2]))

fluxout=ecplotmlr(xvar1=annnet[,i],xvar2=anntrop[,i],xvar3=NULL,yvar=regfluxes[,i],obsval1=obsannnet[,i],obsval2=obsanntrop[,i],obsval3=NULL,nsamp2d=nsample2d,nsampmlr=nsamplemlr,modnames=modnames,expnames=expnames,filenames=c(paste(modeldir,'/NET_Flux_vs_NET-T_Conc_latdiv_MLR_fig4c.png',sep=''),paste(modeldir,'/NET_Flux_vs_NET-T_Conc_latdiv_MLR_figS4c.png',sep='')),plotsel=c(3,5),fitexpmeans=F,inclmodleg=F,inclexpleg=F,ecmain=expression('c. Northern Extratropics'),xlb=expression(paste('Concentration-based Flux Prediction (PgC ',yr^-1,')')),ylb='',xlm=range(regfluxes[,i]),ylm=range(regfluxes[,i]) ) # ,ylb=expression(paste('v10 MIP Total Flux (PgC ',yr^-1,')')) )
remnetlatdiv=c(mean(regfluxes[,i]),sd(regfluxes[,i]),mean(regfluxes[!grepl('OG',expnames),i]),sd(regfluxes[!grepl('OG',expnames),i]),fluxout[1],fluxout[2],fit$coef[1],fit$coef[2],fit$coef[3],NA)
system(paste('cp ',modeldir,'/NET_Flux_vs_NET-T_Conc_latdiv_MLR_fig4c.eps ',modeldir,'/fig4c.eps',sep=''))
system(paste('cp ',modeldir,'/NET_Flux_vs_NET-T_Conc_latdiv_MLR_figS4c.eps ',modeldir,'/figS4c.eps',sep=''))

system(paste('montage ',modeldir,'/SET_Flux_vs_SET-T_Conc_latdiv_MLR_fig4a.png ',modeldir,'/T_Flux_vs_T-ET_Conc_latdiv_MLR_fig4b.png ',modeldir,'/NET_Flux_vs_NET-T_Conc_latdiv_MLR_fig4c.png -geometry +1+3 ',modeldir,'/fig4_ecfcs.png',sep=''))
system(paste('montage ',modeldir,'/SET_Flux_vs_SET-T_Conc_latdiv_MLR_figS4a.png ',modeldir,'/T_Flux_vs_T-ET_Conc_latdiv_MLR_figS4b.png ',modeldir,'/NET_Flux_vs_NET-T_Conc_latdiv_MLR_figS4c.png -geometry +1+3 ',modeldir,'/figS4_ecfcs_bymod.png',sep=''))

# calc simple diff
obsdifmet11=NULL
for(j in c(1:length(sdsel))){
        obsdifmet11=c(obsdifmet11,obsannnet[variables==sdsel[j],i]-obsanntrop[variables==sdsel[j],i])
}
obsdifmet11=c(obsdifmet11,0) # add 0 for NOAA-NOAA
names(obsdifmet11)=c(sdsel,'CO2_NOAA_m_CO2_NOAA')
obsdifremppm=cbind(obsdifremppm,round(obsdifmet11,3))


print('SET-Trop')
print('2-Var MLR')
i=4
fit=lm(regfluxes[,i] ~ annset[,i] + anntrop[,i])
print(paste(labels2[i],round(summary(fit)$adj.r.squared,3),round((1-sd(fit$resid)/sd(regfluxes[,i]))*100,2),round((1-sd(fit$resid)/sd(regfluxes[!grepl('OG',expnames),i]))*100,2)))
obsdifmet=NULL
for(j in c(1:length(sdsel))){
        obsdifmet=c(obsdifmet,(obsannset[variables==sdsel[j],i]*fit$coef[2]+obsanntrop[variables==sdsel[j],i]*fit$coef[3])/mean(abs(fit$coef[2:3])))
}
obsdifmet=c(obsdifmet,0) # add 0 for NOAA-NOAA
obsdifsd=sd(obsdifmet,na.rm=T)
print(obsdifmet)
print(obsdifsd)

# since not using these numbers, need much smaller ensemble
ecplotmlr(xvar1=annset[,i],xvar2=anntrop[,i],xvar3=NULL,yvar=regfluxes[,i],obsval1=obsannset[,i],obsval2=obsanntrop[,i],obsval3=NULL,nsamp2d=nsample2d,nsampmlr=nsamplemlrfast,modnames=modnames,expnames=expnames,filenames=c(paste(modeldir,'/SET-T_Flux_vs_SET-T_Conc_latdiv_MLR.png',sep=''),paste(modeldir,'/SET-T_Flux_vs_SET-T_Conc_latdiv_MLR_bymod.png',sep='')),plotsel=c(3,5),fitexpmeans=F,inclmodleg=F,inclexpleg=F,ecmain=expression('Southern Extratropics - Tropics'),xlb=expression(paste('Concentration-based Flux Prediction (PgC ',yr^-1,')')),ylb='',xlm=range(regfluxes[,i]),ylm=range(regfluxes[,i]) ) 


print('Trop-(SET+NET)')
print('2-Var MLR')
i=5
fit=lm(regfluxes[,i] ~ anntrop[,i] + annet[,i])
print(paste(labels2[i],round(summary(fit)$adj.r.squared,3),round((1-sd(fit$resid)/sd(regfluxes[,i]))*100,2),round((1-sd(fit$resid)/sd(regfluxes[!grepl('OG',expnames),i]))*100,2)))
obsdifmet=NULL
for(j in c(1:length(sdsel))){
        obsdifmet=c(obsdifmet,(obsanntrop[variables==sdsel[j],i]*fit$coef[2]+obsannet[variables==sdsel[j],i]*fit$coef[3])/mean(abs(fit$coef[2:3])))
}
obsdifmet=c(obsdifmet,0) # add 0 for NOAA-NOAA
obsdifsd=sd(obsdifmet,na.rm=T)
print(obsdifmet)
print(obsdifsd)

fluxout=ecplotmlr(xvar1=annet[,i],xvar2=anntrop[,i],xvar3=NULL,yvar=regfluxes[,i],obsval1=obsannet[,i],obsval2=obsanntrop[,i],obsval3=NULL,nsamp2d=nsample2d,nsampmlr=nsamplemlrfast,modnames=modnames,expnames=expnames,filenames=c(paste(modeldir,'/T-ET_Flux_vs_T-ET_Conc_latdiv_MLR.png',sep=''),paste(modeldir,'/T-ET_Flux_vs_T-ET_Conc_latdiv_MLR_bymod.png',sep='')),plotsel=c(3,5),fitexpmeans=F,inclmodleg=F,inclexpleg=F,ecmain=expression('Tropics - Extratropics'),xlb=expression(paste('Concentration-based Flux Prediction (PgC ',yr^-1,')')),ylb='',xlm=range(regfluxes[,i]),ylm=range(regfluxes[,i]) ) 
remetlatdiv=c(mean(regfluxes[,i]),sd(regfluxes[,i]),mean(regfluxes[!grepl('OG',expnames),i]),sd(regfluxes[!grepl('OG',expnames),i]),fluxout[1],fluxout[2],fit$coef[1],fit$coef[2],fit$coef[3],NA)


print('NET-Trop')
print('2-Var MLR')
i=6
fit=lm(regfluxes[,i] ~ annnet[,i] + anntrop[,i])
print(paste(labels2[i],round(summary(fit)$adj.r.squared,3),round((1-sd(fit$resid)/sd(regfluxes[,i]))*100,2),round((1-sd(fit$resid)/sd(regfluxes[!grepl('OG',expnames),i]))*100,2)))
obsdifmet=NULL
for(j in c(1:length(sdsel))){
        obsdifmet=c(obsdifmet,(obsannnet[variables==sdsel[j],i]*fit$coef[2]+obsanntrop[variables==sdsel[j],i]*fit$coef[3])/mean(abs(fit$coef[2:3])))
}
obsdifmet=c(obsdifmet,0) # add 0 for NOAA-NOAA
obsdifsd=sd(obsdifmet,na.rm=T)
print(obsdifmet)
print(obsdifsd)

ecplotmlr(xvar1=annnet[,i],xvar2=anntrop[,i],xvar3=NULL,yvar=regfluxes[,i],obsval1=obsannnet[,i],obsval2=obsanntrop[,i],obsval3=NULL,nsamp2d=nsample2d,nsampmlr=nsamplemlrfast,modnames=modnames,expnames=expnames,filenames=c(paste(modeldir,'/NET-T_Flux_vs_NET-T_Conc_latdiv_MLR.png',sep=''),paste(modeldir,'/NET-T_Flux_vs_NET-T_Conc_latdiv_MLR_bymod.png',sep='')),plotsel=c(3,5),fitexpmeans=F,inclmodleg=F,inclexpleg=F,ecmain=expression('Northern Extratropics - Tropics'),xlb=expression(paste('Concentration-based Flux Prediction (PgC ',yr^-1,')')),ylb='',xlm=range(regfluxes[,i]),ylm=range(regfluxes[,i]) )


# Print output

outputlatdiv=cbind(remsetlatdiv,remtlatdiv,remnetlatdiv,remetlatdiv)
outputlatdiv=rbind(outputlatdiv,(1-outputlatdiv[6,]/outputlatdiv[2,])*100,(1-outputlatdiv[6,]/outputlatdiv[4,])*100) 
rownames(outputlatdiv)=c('Prior Flux','Prior Flux SD','Prior Flux no OG','Prior Flux SD no OG','Flux result','Flux result SD','Fit Coef 1','Fit Coef 2','Fit Coef 3','Fit Coef 4','Uncred','Uncred no OG') # ,'Unc % All Exp','Unc % no OG','Unc % EC')
print(round(outputlatdiv,2))
print(latdivs)

fitcoef=outputlatdiv[7:8,]
alpha=2*outputlatdiv[7,]/(abs(outputlatdiv[7,])+abs(outputlatdiv[8,]))
beta=-2*outputlatdiv[8,]/(abs(outputlatdiv[7,])+abs(outputlatdiv[8,]))
print(rbind(fitcoef,alpha,beta))
print('old weighted gradients')
predicted_gradients=c(obsconc[1]*alpha[1]-obsconc[2]*beta[1],obsconc[2]*alpha[2]-obsconc[4]*beta[2],obsconc[3]*alpha[3]-obsconc[2]*beta[3])
print(predicted_gradients)
print('conc-mlo')
print(obsconc)
print('SET-T, T-ET, NET-T')
print(c(obsconc[1]-obsconc[2],obsconc[2]-obsconc[4],obsconc[3]-obsconc[2]))

print(obsdifremppm)
print(apply(obsdifremppm,2,sd))
print(obsdifrempgc)
print(apply(obsdifrempgc,2,sd))

# print results in Latex format for tables

table2remcoef=round(table2remcoef,2)
table2remse=round(table2remse,2)
print('Table S2 (coef)')
print(paste('$F_{NET}$ \newline $>$ 20$^\\circ$ N & {} & ',table2remcoef[3,2],' $\\pm$ ',table2remse[3,2],' & ',table2remcoef[3,3],' $\\pm$ ',table2remse[3,3],' & ',table2remcoef[3,4],' $\\pm$ ',table2remse[3,4],' \\',sep=''))
print(paste('$F_T$ \newline 20$^\\circ$ S - 20$^\\circ$ N & ',table2remcoef[2,1],' $\\pm$ ',table2remse[2,1],' & ',table2remcoef[2,2],' $\\pm$ ',table2remse[2,2],' & ',table2remcoef[2,3],' $\\pm$ ',table2remse[2,3],' & ',table2remcoef[2,4],' $\\pm$ ',table2remse[2,4],' \\',sep=''))
print(paste('$F_{SET}$ \newline $>$ 20$^\\circ$ S & ',table2remcoef[1,1],' $\\pm$ ',table2remse[1,1],' & ',table2remcoef[1,2],' $\\pm$ ',table2remse[1,2],' & {} & ',table2remcoef[1,4],' $\\pm$ ',table2remse[1,4],' \\',sep=''))
print('Remove extra slash from in front of circ and pm')

obsdifmetglb=NULL
for(j in c(1:length(sdsel))){
        obsdifmetglb=c(obsdifmetglb,obsannglb[variables==sdsel[j],i])
}
obsdifmetglb=c(obsdifmetglb,0)
obsdifremppm=cbind(obsdifremppm,obsdifmetglb)

print('Table 1 (results)')
print(paste('v10 MIP total flux \\newline full ensemble & ',round(outputlatdiv[1,1],2),' $\\pm$ ',round(outputlatdiv[2,1],2),' & ',round(outputlatdiv[1,2],2),' $\\pm$ ',round(outputlatdiv[2,2],2),' & ',round(outputlatdiv[1,3],2),' $\\pm$ ',round(outputlatdiv[2,3],2),' \\',sep=''))
print(paste('v10 MIP total flux \\newline no OG subset & ',round(outputlatdiv[3,1],2),' $\\pm$ ',round(outputlatdiv[4,1],2),' & ',round(outputlatdiv[3,2],2),' $\\pm$ ',round(outputlatdiv[4,2],2),' & ',round(outputlatdiv[3,3],2),' $\\pm$ ',round(outputlatdiv[4,3],2),' \\',sep=''))
print(paste('\textbf{ATom + v10 MIP \\newline  ECFC total flux} & \textbf{',round(outputlatdiv[5,1],2),' $\\pm$ ',round(outputlatdiv[6,1],2),'} & \textbf{',round(outputlatdiv[5,2],2),' $\\pm$ ',round(outputlatdiv[6,2],2),'} & \textbf{',round(outputlatdiv[5,3],2),' $\\pm$ ',round(outputlatdiv[6,3],2),'} \\',sep=''))
print(paste('Uncertainty reduction \\newline from full ensemble (\\%) & ',round(outputlatdiv[11,1],1),' & ',round(outputlatdiv[11,2],1),' & ',round(outputlatdiv[11,3],1),' \\',sep=''))
print(paste('Uncertainty reduction \\newline from no OG subset (\\%)  & ',round(outputlatdiv[12,1],1),' & ',round(outputlatdiv[12,2],1),' & ',round(outputlatdiv[12,3],1),' \\',sep=''))
print('remove slashes from \\newline, \\pm, and \\%, and add zeros')

## to copy to foss_ocean.r:
print(paste('netecfc=',round(outputlatdiv[5,3],3),'; uncnetecfc=',round(outputlatdiv[6,3],3),sep=''))
print(paste('trpecfc=',round(outputlatdiv[5,2],3),'; unctrpecfc=',round(outputlatdiv[6,2],3),sep=''))
print(paste('setecfc=',round(outputlatdiv[5,1],3),'; uncsetecfc=',round(outputlatdiv[6,1],3),sep=''))


print('Table S3 (ppm)')
print(paste('Observed Gradient & ',paste(round(c(obsconc[1]-obsconc[2],obsconc[2]-obsconc[4],obsconc[3]-obsconc[2],obsconc[5]),2),collapse=' & '),' \\',sep='')) 
print(paste('Harvard QCLS & ',paste(round(obsdifremppm[2,],2),collapse=' & '),' \\',sep='')) 
print(paste('NSF NCAR AO2 & ',paste(round(obsdifremppm[1,],2),collapse=' & '),' \\',sep='')) 
print(paste('NSF NCAR / Scripps Medusa & ',paste(round(obsdifremppm[3,],2),collapse=' & '),' \\',sep='')) 
print(paste('NOAA PFP & ',paste(round(obsdifremppm[4,],2),collapse=' & '),' \\',sep='')) 
print('NOAA Picarro & 0.00 & 0.00 & 0.00 & 0.00 \\')
print(paste('1-$\\sigma$ & ',paste(round(apply(obsdifremppm,2,sd),2),collapse=' & '),' \\',sep='')) 
print('Remove extra slash from in front of sigma')

print('Table S4 (PgC)')
print(paste('Harvard QCLS & ',paste(round(obsdifrempgc[2,],2),collapse=' & '),' \\',sep='')) 
print(paste('NSF NCAR AO2 & ',paste(round(obsdifrempgc[1,],2),collapse=' & '),' \\',sep='')) 
print(paste('NSF NCAR / Scripps Medusa & ',paste(round(obsdifrempgc[3,],2),collapse=' & '),' \\',sep='')) 
print(paste('NOAA PFP & ',paste(round(obsdifrempgc[4,],2),collapse=' & '),' \\',sep='')) 
print('NOAA Picarro & 0.00 & 0.00 & 0.00 \\')
print(paste('1-$\\sigma$ & ',paste(round(apply(obsdifrempgc,2,sd),2),collapse=' & '),' \\',sep='')) 
print('Remove extra slash from in front of sigma')

print('T-ET results')
# these are for T-ET conc vs T-ET flux, whereas others are not diff flux (e.g. remnetlatdiv is just T-ET conc vs NET flux)
print(paste('v10 MIP total flux \\newline full ensemble & ',round(outputlatdiv[1,4],2),' $\\pm$ ',round(outputlatdiv[2,4],2),sep=''))
print(paste('v10 MIP total flux \\newline no OG subset & ',round(outputlatdiv[3,4],2),' $\\pm$ ',round(outputlatdiv[4,4],2),sep=''))
print(paste('\textbf{ATom + v10 MIP \\newline  ECFC total flux} & \textbf{',round(outputlatdiv[5,4],2),' $\\pm$ ',round(outputlatdiv[6,4],2),sep=''))
print(paste('Uncertainty reduction \\newline from full ensemble (\\%) & ',round(outputlatdiv[11,4],1),sep=''))
print(paste('Uncertainty reduction \\newline from no OG subset (\\%)  & ',round(outputlatdiv[12,4],4),sep=''))
print('remove slashes from \\newline, \\pm, and \\%, and add zeros')
