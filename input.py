########################################
# create input data file
# input.py
########################################

import numpy as np

#ac=30. 
ms=1.
md=0.1

# Cambiar A = 0.0
# A = 0.1 f = 0.3
# A = 0.7 f = 1

Apert=0.7 
Fpert=1.

gama=1.
cmigI=0.1

#metal=2.069580554962158E-002
metal=0.0
tgas=6859750.19

Tfin=3.e7
Verb=False

acvals=np.linspace(20,70,300)

f = open('parameters.in', 'w')
f.write("md,ac,ms,metal,tgas,gama,cmigI,Apert,Fpert,Tfin,Verb\n") 
for ac in acvals:
    amin=ac*((7/4.-gama)/(2.-gama))**(1./(2-gama))
    sigmag=(2.-gama)*md*2.e33/(2.*np.pi*(ac*1.5e13)**2.)
    qc=1.24e5*(amin**(gama-7./4.))*(ac**(-gama))*(ms)*np.exp((amin/ac)**(2.-gama))/sigmag
    if(qc>1):
        f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(md,ac,ms,metal,tgas,gama,cmigI,Apert,Fpert,Tfin,Verb)) 
f.close()
