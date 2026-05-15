source('harm_fit.r') # harmonic fitting code

source('harm_plot.r') # harmonic fitting code

firstrun=T # can set to F after running once to avoid re-reading ascii 10-sec mergefiles
if(firstrun){
	source('read_airborne_obs.r')
	readac() # outdir='OBS', inclmed=T, and inclpfp=T in defaults
}

# Set default parameters for all calls to harmf and harmp
obsoutdir='OBS/BINFIT'; modoutdir='V10MIP'

hijack <- function (FUN, ...) {
    .FUN <- FUN
    args <- list(...)
    invisible(lapply(seq_along(args), function(i) {
        formals(.FUN)[[names(args)[i]]] <<- args[[i]]
    }))
    .FUN
}

.harmf=hijack(harmf,latspan=10,prsspan=100,flag='_10x100',latrange=c(-70,80),prsrange=c(300,1000),aggby='campflt',adjMLO2X2007=T,units='ppm',binxaxcex=0.75,binyaxcex=0.75,binprscex=0.7,binlatcex=0.6)
## for OCO-2 comparison, use nharm=1, useallrf=F (only over ocean, excludes A1 RFs over NAm), inclmed=F
# adjusting MLO back to X2007 as ATom obs (incl local AO2 and MED used here) and v10MIP both on X2007

.harmp=hijack(harmp,flag='_10x100',xlim1=c(-70,80),ylim1=c(1000,300),units='ppm',pngh=500)


## First process obs 

print('Fitting and plotting observations')

# THETA
print('THETA')
.harmf(species='THETA',units='K',detrend=F,latrange=c(-70,80),outdir=obsoutdir)
.harmp(species='THETA',specname='Theta',units='K',detrend=F,xlim1=c(-70,80),inputdir=obsoutdir,netdir=obsoutdir,zlimann=c(-4.6,2.6),zlimseas=c(0,20),zlimmon=c(-11,7),maintf=T)

# CO2
variables=c('CO2_OP','CO2_AO2_m_CO2_NOAA','CO2_QCLS_m_CO2_NOAA','CO2_MED_m_CO2_NOAA','CO2_PFP_m_CO2_NOAA','CO2_NOAA') # added CO2_NOAA to allow comparison
# ObsPack is sampling CO2_NOAA except first 2 flights of ATom-1 - CO2.X is similar but includes QCLS in cal gaps

for(var in variables){

	print(var)

	if(grepl('CO2_MED',var)){ 
		.harmf(species=var,detrend=!grepl('_m_',var),outdir=obsoutdir,tomfilenamein='RData.ATom.Mer.all_MED_trimmed')
	} else if(grepl('CO2_PFP',var)){ 
		.harmf(species=var,detrend=!grepl('_m_',var),outdir=obsoutdir,tomfilenamein='RData.ATom.Mer.all_PFP_trimmed')
	} else {
		.harmf(species=var,detrend=!grepl('_m_',var),outdir=obsoutdir)
	}

	.harmp(species=var,specname=var,detrend=!grepl('_m_',var),inputdir=obsoutdir,netdir=obsoutdir,maintxt='Observations')

}


## Then process OCO2-MIP models

# V10MIP ATom ObsPack data already merged in with flight data by model_aircraft_merge.r,
# which writes files ATom_Merged_Model_Output.txt in MOD/EXP subdirectories

print('Fitting and plotting OCO2-MIP models')

# Options:
models=c('AMES','BAKER','CAMS','CMS-FLUX','COLA','CT','OU','TM5-4DVAR','UT','WOMBAT')
experiments=c('IS','OG','LNLG','LNLGOGIS','LNLGIS')

for(exp in experiments){
        print(exp)

        for(mod in models){

		print(mod)

		.harmf(species='CO2_OP',tomfilenamein='ATom_Merged_Model_Output.txt',oco2mipoutput=T,indir=paste(modoutdir,'/',mod,'/',exp,sep=''),outdir=paste(modoutdir,'/',mod,'/',exp,sep=''))

		.harmp(species='CO2_OP',specname='CO2',inputdir=paste(modoutdir,'/',mod,'/',exp,sep=''),
		       netdir=paste(modoutdir,'/',mod,'/',exp,sep=''),maintxt=paste(mod,exp)) 

	}

}
