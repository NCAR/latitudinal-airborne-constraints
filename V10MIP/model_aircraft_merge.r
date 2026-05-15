# Program to read in OCO2 V10MIP model ATom obspack files and merge onto selected 10-sec aircraft variables

## also calculates cross-model variability
## also reads in observed CO2 and calculates median model offsets
## have to run ../OBS/aircraft_obspack_merge.r first

## V10MIP was all on WMO X2007 scale


library('ncdf4')

## Options:

modobspackdir='../LOCALDATA/V10MIP'
op='obspack_GLOBALVIEWplus_v6.1'
obspackdir='../LOCALDATA' # still on X2007

models=c('AMES','BAKER','CAMS','CMS-FLUX','COLA','CT','OU','TM5-4DVAR','UT','WOMBAT') # not using CSU, JHU, NIES, or LoFi
experiments=c('IS','OG','LNLG','LNLGOGIS','LNLGIS','unopt')

## Process concentrations:

# Open single netcdf file for all models/experiments ATom results:
print('reading in ATom model output ObsPack format')
modtomnc=nc_open(paste(modobspackdir,'/co2_tom_aircraft-insitu_1_allvalid.nc',sep=''))
opncsubs=ncvar_get(modtomnc,'submissions')
#> opncsubs
#     "CT"        "CAMS"      "UT"        "Baker"     "WOMBAT"    "CSU" #     "TM5-4DVAR" "OU"        "Ames"      "COLA"      "CMS-Flux"  "LoFI"       "M2CC"      "JHU"       "CSU-NEE"   "NIES"
opncexps=ncvar_get(modtomnc,'experiments')
#> opncexps
# "IS"       "LNLG"     "OG"       "LNLGIS"   "LNLGOGIS" "unopt"
opncids=ncvar_get(modtomnc,'obspack_id') # not using
opncval=ncvar_get(modtomnc,'simulated_values')
## V10 missing CAMS unopt, CSU all, JHU all, CSU-NEE all, NIES IS and LNLG (and NIES for others missing 90%)

# Open observation ObsPack too
print('reading in ATom obs data ObsPack')
obstomnc=nc_open(paste(obspackdir,'/co2_tom_aircraft-insitu_1_allvalid.nc',sep=''))

# Subselect ATom model output position data and observed CO2
modtompos=data.frame(cbind(t(ncvar_get(modtomnc,'time_components')),ncvar_get(modtomnc,'altitude'),ncvar_get(modtomnc,'latitude'),ncvar_get(modtomnc,'longitude')),ncvar_get(modtomnc,'obspack_id')) ; colnames(modtompos)=c('year','mon','day','hour','min','sec','alt','lat','lon','obspack_id')
nc_close(modtomnc)

obstompos=data.frame(cbind(t(ncvar_get(obstomnc,'time_components')),ncvar_get(obstomnc,'altitude'),ncvar_get(obstomnc,'latitude'),ncvar_get(obstomnc,'longitude')),ncvar_get(obstomnc,'obspack_id')) ; colnames(obstompos)=c('year','mon','day','hour','min','sec','alt','lat','lon','obspack_id')
obsco2=ncvar_get(obstomnc,'value')*1E6 # X2007 for GV+ 6.1 and earlier
nc_close(obstomnc)

# Read in merged ObsPack and 10-sec data
mrgtomdat=read.table(paste('../OBS/ATom_',op,'_merge.txt',sep=''),header=T,stringsAsFactors=F)

## there more records in the OpsPack than in the 10-sec merge product, mostly because of test flight data included in the ObsPack

# Merge in ancillary flight data
varstoadd=c('camp','flt','prof','pressure','theta','strat','hdg')
modtompos=cbind(modtompos,matrix(NA,nrow(modtompos),length(varstoadd))); colnames(modtompos)=c('year','mon','day','hour','min','sec','alt','lat','lon','obspack_id',varstoadd)

modtompos[is.element(modtompos$obspack_id,mrgtomdat$obspack_id),varstoadd]=mrgtomdat[match(modtompos$obspack_id[is.element(modtompos$obspack_id,mrgtomdat$obspack_id)],mrgtomdat$obspack_id),varstoadd] # not replacing lat, lon, alt
print(paste('Models missing ATom merged camp, flt, prof, pressure, theta, strat, hdg =',sum(is.na(modtompos$camp)),' out of',nrow(modtompos)))

obstompos[is.element(obstompos$obspack_id,mrgtomdat$obspack_id),varstoadd]=mrgtomdat[match(obstompos$obspack_id[is.element(obstompos$obspack_id,mrgtomdat$obspack_id)],mrgtomdat$obspack_id),varstoadd] # not replacing lat, lon, alt
print(paste('Obs missing ATom merged camp, flt, prof, pressure, theta, strat =',sum(is.na(obstompos$camp)),' out of',nrow(obstompos)))

tomdatdt=as.POSIXlt(ISOdatetime(modtompos$year,modtompos$mon,modtompos$day,modtompos$hour,modtompos$min,modtompos$sec,tz='UTC'))
obstomdatdt=as.POSIXlt(ISOdatetime(obstompos$year,obstompos$mon,obstompos$day,obstompos$hour,obstompos$min,obstompos$sec,tz='UTC'))

print('Before filtering model aircraft files')
print(dim(modtompos))
print(dim(opncval))
print(length(tomdatdt))
print(dim(obstompos))
print(length(obstomdatdt))

modtomposunfilt=modtompos
opncvalunfilt=opncval
tomdatdtunfilt=tomdatdt
obsco2unfilt=obsco2

# Calculate cross-model SDs by experiment
modsd=apply(opncvalunfilt,c(1,3),sd,na.rm=T)

# Calculate absolute median mod-obs differences by experiment
medmodminobs=apply(opncvalunfilt,c(1,3),median,na.rm=T)
for(i in c(1:6)){ # 6 is prior
	medmodminobs[i,]=abs(medmodminobs[i,]-obsco2unfilt)
}


# Loop on experiment and model:

# Write out mod/exp-specific files and one composite file
mtpnames=colnames(modtompos)
alltomdat=cbind(modtompos,obsco2,modsd[4,],medmodminobs[4,])
colnames(alltomdat)=c(mtpnames,'obs','modsd','medmodminobs')

for(exp in c(1:length(experiments))){ # All
	exper=experiments[exp]
	print(exper)

	for(mod in c(1:length(models))){ # All

		model=models[mod] # dirs are all caps
		modelname=model # actual name used in netcdf file
		if(modelname=='CMS-FLUX') modelname='CMS-Flux'; if(modelname=='AMES') modelname='Ames'; if(modelname=='BAKER') modelname='Baker'; if(modelname=='WEIR') modelname='LoFI'
		if(exper!='LNLGOGIS'|model!='CSU'){

			if(modelname!='LoFI'|exper=='IS'){

				print(model)
				
				tomdat=cbind(modtompos,opncval[which(opncexps==exper),which(opncsubs==modelname),])
				colnames(tomdat)[ncol(tomdat)]='co2'
				alltomdat=cbind(alltomdat,opncval[which(opncexps==exper),which(opncsubs==modelname),])
				colnames(alltomdat)[ncol(alltomdat)]=paste(model,'_',exper,sep='')

				if(!file.exists(model)) system(paste('mkdir',model))
				if(!file.exists(paste(model,'/',exper,sep=''))) system(paste('mkdir ',model,'/',exper,sep=''))
				write(colnames(tomdat),paste(model,'/',exper,'/ATom_Merged_Model_Output.txt',sep=''),ncol=ncol(tomdat))
				write(t(tomdat),paste(model,'/',exper,'/ATom_Merged_Model_Output.txt',sep=''),ncol=ncol(tomdat),append=T)

			} # excluding LoFI for all but IS

		} # excluding CSU/LNLGOGIS

	} # loop on model

} # loop on experiment


write(colnames(alltomdat),'All_ATom_Merged_Model_Output.txt',ncol=ncol(alltomdat))
write(t(alltomdat),'All_ATom_Merged_Model_Output.txt',ncol=ncol(alltomdat),append=T)
