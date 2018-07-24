# correr como python hacer_animacion.py archivo.txt
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
filename=sys.argv[1]
command = 'awk \'{{print $1}}\' {} | uniq -c > indices.txt'.format(filename)
os.system(command)
ind=np.loadtxt('indices.txt',skiprows=1)
rawdata=np.loadtxt(filename,skiprows=1)
ncols=6
myrange=len(ind)
figure_list = ''
ii=0
trange=int(ind[:myrange,0].sum())
radmax=rawdata[:trange,5].max()
radmax2=rawdata[:,5].max()
#eradcm=6.37e8
eradcm=1
tming=5e6
mming=15
figure_list2 = ''
giant=False
if rawdata[:,4].max() > mming:
    giant=True
for j in np.arange(myrange):
    plt.xlim(0,rawdata[:trange,2].max()+0.2)
# cambiar aqui para tener eje y segun el planeta mas grande o fijo
    plt.ylim(0,rawdata[:trange,4].max()+0.2)
    plt.ylim(0,3)
    data=np.zeros((int(ind[j,0]),ncols))
    for i in range(int(ind[j,0])):
        for k in range(ncols):
            data[i,k]=rawdata[ii,k]
        ii=ii+1
    plt.scatter(data[:,2],data[:,4],s=100*data[:,5],c=data[:,5]/eradcm,cmap='cubehelix',vmin=0,vmax=radmax/eradcm)   
    plt.colorbar(label=r'Planet radius ($R_\oplus$)');
    plt.xlabel('Semi-major axis (AU)')
    plt.ylabel(r'Planet mass ($M_\oplus$)')
    plt.title('Time = %.1e yr' % data[i,1])
    figure_name='plot%06d.png' % j
    plt.savefig(figure_name)
    figure_list = figure_list + ' {} '.format(figure_name)
    plt.close() 
    if giant and (data[0,1] > tming):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax2=ax.scatter(data[:,2],data[:,4],s=100*data[:,5],c=data[:,5]/eradcm,cmap='cubehelix',vmin=0,vmax=radmax2/eradcm)   
        fig.colorbar(ax2,label=r'Planet radius ($R_\oplus$)')
        ax.set_xlim(0,rawdata[:trange,2].max()/2.)
        ax.set_ylim(0,rawdata[:trange,4].max()+2)
        ax.set_xlabel('Semi-major axis (AU)')
        ax.set_ylabel(r'Planet mass ($M_\oplus$)')
        ax.set_title('Time = %.1e yr' % data[i,1])
        figure_name2='gplot%06d.png' % j
        fig.savefig(figure_name2)
        figure_list2 = figure_list2 + ' {} '.format(figure_name2)
        plt.close()

for i in range(15):    
    figure_list = figure_list+ ' {} '.format(figure_name)  
if giant and (data[0,1] > tming):    
	for i in range(15):    
	    figure_list2 = figure_list2+ ' {} '.format(figure_name2)  
	command2 = 'convert -delay 10 -loop 0 {}'.format(figure_list2)+' ganimation_{}.gif'.format(filename[:-4])
	os.system(command2)

command = 'convert -delay 10 -loop 0 {}'.format(figure_list)+' animation_{}.gif'.format(filename[:-4])
os.system(command)
command = 'rm *.png'
os.system(command)
command = 'rm indices.txt'
os.system(command)
