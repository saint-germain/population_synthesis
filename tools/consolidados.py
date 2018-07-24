import numpy as np
import matplotlib.pyplot as plt
import os

fig = plt.figure(figsize=(24,20))
fig2 = plt.figure(figsize=(30,20))
fig3 = plt.figure(figsize=(24,20))

fig4 = plt.figure(figsize=(24,20))
fig5 = plt.figure(figsize=(30,20))

Alist=np.array([0.1,0.3,0.5,0.7])
flist=np.array([0.3,0.7,1.,3.])
XX, YY = np.meshgrid(Alist, flist)
XX=XX.T
YY=YY.T
filename="Original50mopt5.txt"
command = 'awk \'{{print $1}}\' {} | uniq -c > indices.txt'.format(filename)
os.system(command)
ind=np.loadtxt('indices.txt',skiprows=1)
rawdata=np.loadtxt(filename,skiprows=1)
ncols=6
myrange=len(ind)
#myrange=200
figure_list = ''
ii=0
trange=int(ind[:myrange,0].sum())
radmax=rawdata[:trange,5].max()
eradcm=6.37e8
armass=[]
agmass=[]
acmass=[]
anplan=[]
atime=[]
for j in np.arange(myrange):
    data=np.zeros((int(ind[j,0]),ncols))
    for i in range(int(ind[j,0])):
        for k in range(ncols):
            data[i,k]=rawdata[ii,k]
        ii=ii+1
    rmass=data[:,4].sum()      #rocky mass
    gmass=data[:,3].sum()    #gas mass
    cmass=(data[:,2]*(data[:,4]+data[:,3])).sum()/(rmass+gmass)     # center of mass
    nplan=len(data)*1.
    armass.append(rmass)
    agmass.append(gmass)
    acmass.append(cmass)
    anplan.append(nplan)
    atime.append(data[0,1])
armass=np.asarray(armass)
agmass=np.asarray(agmass)
acmass=np.asarray(acmass)
anplan=np.asarray(anplan)
atime=np.asarray(atime)
armass2=armass
agmass2=agmass
acmass2=acmass
anplan2=anplan
atime2=atime
for kl in range(1,17):
    filename = "pert50M%02da.txt" % kl
#    filename="perturbacion01.txt"
    command = 'awk \'{{print $1}}\' {} | uniq -c > indices.txt'.format(filename)
    os.system(command)
    ind=np.loadtxt('indices.txt',skiprows=1)
    rawdata=np.loadtxt(filename,skiprows=1)
    ncols=6
    myrange=len(ind)
    #myrange=200
    figure_list = ''
    ii=0
    trange=int(ind[:myrange,0].sum())
    radmax=rawdata[:trange,5].max()
    eradcm=6.37e8
    armass=[]
    agmass=[]
    acmass=[]
    anplan=[]
    atime=[]
    for j in np.arange(myrange):
        data=np.zeros((int(ind[j,0]),ncols))
        for i in range(int(ind[j,0])):
            for k in range(ncols):
                data[i,k]=rawdata[ii,k]
            ii=ii+1
        rmass=data[:,4].sum()      #rocky mass
        gmass=data[:,3].sum()    #gas mass
        cmass=(data[:,2]*(data[:,4]+data[:,3])).sum()/(rmass+gmass)     # center of mass
        nplan=len(data)*1.
        armass.append(rmass)
        agmass.append(gmass)
        acmass.append(cmass)
        anplan.append(nplan)
        atime.append(data[0,1])
    armass=np.asarray(armass)
    agmass=np.asarray(agmass)
    acmass=np.asarray(acmass)
    anplan=np.asarray(anplan)
    atime=np.asarray(atime)
    
    intmass=armass+agmass
    intmass2=armass2+agmass2
    ax = fig.add_subplot(4,4,kl)
    ax.plot(atime,anplan,label="With perturbation") #this one side to side with the next two-in one
    ax.plot(atime2,anplan2,'r--',label="No perturbation")
    ax.set_xlabel(r'$t$ (yr)' )
    ax.set_ylabel(r'$N_p$' )
    ax.set_title(r'$A\ =\ %.1f$, $f\ =\ %.1f$'%(np.reshape(XX,16)[kl-1],np.reshape(YY,16)[kl-1]))
    figure_name1="r_nplanet.png"
    ax.legend()
    fig.savefig(figure_name1)

    ax2 = fig2.add_subplot(4,4,kl)
    ax2.plot(atime,(intmass)/anplan,lw=2,label="With perturbation")
    ax2.set_xlabel(r'$t$ (yr)' )
    ax2.set_ylabel(r'$M_p/N_p\ (M_\oplus)$' )
    ax2.set_title(r'$A\ =\ %.1f$, $f\ =\ %.1f$'%(np.reshape(XX,16)[kl-1],np.reshape(YY,16)[kl-1]))
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
    axx2 = ax2.twinx()
    axx2.plot(atime,intmass/anplan,c='r',label="With perturbation")    
    axx2.plot(atime2,intmass2/anplan2,'r--',label="No perturbation")
    axx2.set_ylim(0,(intmass2/anplan2).max()+1)
    for tl in axx2.get_yticklabels():
        tl.set_color('r')
    figure_name2="r_movernplanet.png"
    ax2.legend()
    fig2.savefig(figure_name2)
        
    ax3=fig3.add_subplot(4,4,kl)
    ax3.plot(atime,agmass,label="With perturbation")
    ax3.plot(atime2,agmass2,'r--',label="No perturbation")    
    ax3.set_xlabel(r'$t$ (yr)' )
    ax3.set_ylabel(r'Gas mass ($M_\oplus$)' )
    ax3.set_xlim(0,5e7)
    ax3.set_title(r'$A\ =\ %.1f$, $f\ =\ %.1f$'%(np.reshape(XX,16)[kl-1],np.reshape(YY,16)[kl-1]))
    figure_name3="r_mgas.png"
    ax3.legend()
    fig3.savefig(figure_name3)

    ax4 = fig4.add_subplot(4,4,kl)
    ax4.plot(atime,acmass,label="With perturbation")
    ax4.plot(atime2,acmass2,'r--',label="No perturbation")    
    ax4.set_xlabel(r'$t$ (yr)' )
    ax4.set_ylabel(r'Center of Mass (AU)' )
    ax4.set_xlim(0,5e7)
    ax4.set_title(r'$A\ =\ %.1f$, $f\ =\ %.1f$'%(np.reshape(XX,16)[kl-1],np.reshape(YY,16)[kl-1]))
    figure_name4="r_centerofmass.png"
    ax4.legend()
    fig4.savefig(figure_name4)
    
    ax5 = fig5.add_subplot(4,4,kl)
    ax5.plot(atime,intmass,lw=2,label="With perturbation")
    ax5.set_xlabel(r'$t$ (yr)')
    ax5.set_ylabel(r'$M_p$ $(M_\oplus)$')
    for tl in ax5.get_yticklabels():
        tl.set_color('b')
    axx2 = ax5.twinx()
    axx2.plot(atime,intmass,c='r',label="With perturbation")    
    axx2.plot(atime2,intmass2,'r--',label="No perturbation")
    axx2.set_ylim(0,intmass2.max()+5)
    for tl in axx2.get_yticklabels():
        tl.set_color('r')
    figure_name5="r_mplanet.png"
    ax5.legend()
    fig5.savefig(figure_name5)
