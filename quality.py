# Quick code for evaluating quality of simulations
import numpy as np
import sys
cre_f=sys.argv[1]
tot_f=sys.argv[2]
fin_f=sys.argv[3]
cre=np.loadtxt(cre_f)
tot=np.loadtxt(tot_f)
fin=np.loadtxt(fin_f)
print "*****************************"
print "Planned simulations = "+str(tot.shape[0])
print "Saved to "+tot_f
print "*****************************"
ne=tot[np.logical_not(np.in1d(tot,cre[:,0]))] # not executed
ne_f="notexec.txt"
np.savetxt(ne_f,ne,fmt="%4i")
print "Not executed (interrupted) = "+str(ne.shape[0])
print "Saved to "+ne_f
print "*****************************"
ioe=cre[np.logical_not(np.in1d(cre[:,0],fin))][:,0] # executed but not finished - interrupted or error
ioe_f="intorerr.txt"
np.savetxt(ioe_f,ioe,fmt="%4i")
print "Executed but interrupted (or error) = "+str(ioe.shape[0])
print "Saved to "+ioe_f
print "*****************************"
print "Successful simulations (finished) = "+str(fin.shape[0])
print "Saved to "+fin_f
print "*****************************"
print "Sanity check, PS=NE+IE+FS: "+str(tot.shape[0])+" = "+str(ne.shape[0]+ioe.shape[0]+fin.shape[0])
print "*****************************"
ef=cre[cre[:,1]==0][:,0] # empty files
fs=ef[np.logical_not(np.in1d(ef,ioe))] # failed systems
fs_f="failedsys.txt"
np.savetxt(fs_f,fs,fmt="%4i")
print "Failed systems (no planets formed) = "+str(fs.shape[0])
print "Saved to "+fs_f
print "*****************************"
