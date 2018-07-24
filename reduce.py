import pandas as pd
import numpy as np
import sys
resfile=sys.argv[1]
dff=pd.read_csv(resfile)
header=np.array(["ident, com, nplanets, massbudget, mass efficiency, sigmag0, md, rc, ms"])
np.savetxt("header.txt",header,fmt='%s')
dfter=dff[dff['emepla(i)/emet']<10.]
rdata=np.zeros((len(np.unique(dfter.ident)),9))
kk=0
for i in np.unique(dfter.ident):    
    filter=dfter.ident==i
    com=((dfter[filter]['emepla(i)/emet']*dfter[filter]['a(i)']).sum())/dfter[filter]['emepla(i)/emet'].sum()
    npl=filter.sum()
    mbud=(dfter[filter]['emepla(i)/emet']).sum()
    effm=mbud*3e-6/(dfter.emed[filter].iloc[0])
    sigmag0=dfter.sigmag_0[filter].iloc[0]
    md=dfter.emed[filter].iloc[0]
    rc=dfter.rc[filter].iloc[0]
    ms=dfter.emestar[filter].iloc[0]
    rdata[kk,:]=i,com,npl,mbud,effm,sigmag0,md,rc,ms    
    kk=kk+1
np.savetxt("terrestrial.txt",rdata)
dfgia=dff[dff['emepla(i)/emet']>10.]
rdatag=np.zeros((len(np.unique(dfgia.ident)),9))
kk=0
for i in np.unique(dfgia.ident):    
    filter=dfgia.ident==i
    com=((dfgia[filter]['emepla(i)/emet']*dfgia[filter]['a(i)']).sum())/dfgia[filter]['emepla(i)/emet'].sum()
    npl=filter.sum()
    mbud=(dfgia[filter]['emepla(i)/emet']).sum()
    effm=mbud*3e-6/(dfgia.emed[filter].iloc[0])
    sigmag0=dfgia.sigmag_0[filter].iloc[0]
    md=dfgia.emed[filter].iloc[0]
    rc=dfgia.rc[filter].iloc[0]
    ms=dfgia.emestar[filter].iloc[0]
    rdatag[kk,:]=i,com,npl,mbud,effm,sigmag0,md,rc,ms    
    kk=kk+1
np.savetxt("giant.txt",rdatag)
print("Total systems: "+str(len(np.unique(dff.ident))))
print("Systems with giant planets: "+str(kk-1))
print("%.1f percent of all planets are below 10 earth masses" %(100.*len(dff[dff['emepla(i)/emet']<10])/len(dff)))
