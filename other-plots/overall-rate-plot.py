#!/usr/bin/env python
import sys
import os
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import pdb
import warnings


def portray(center,low,high):
        expo=int(np.floor(np.log10(center)))
        mcenter=center/10**expo
        #dhi=(high-center)/10**expo
        #dlo=(center-low)/10**expo
        mhi=high/10**expo
        mlo=low/10**expo
        s='{0:.2f} \\times 10${{}}^{3:d}$ & {1:.2f}-{2:.2f} \\times 10${{}}^{3:d}$'.format(mcenter,mlo,mhi,expo)
        return s

#format: all-data est., median est., 68% CI lo, 68% CI hi, all in s^-1
computational_rates={'1':(4.28,2.78,0.77,9.07),'2 (AMBER)':(2.06e4,1.9e4,8.82e3,4.05e4),
	'2 (CHARMM)':(3.03e3,2.79e3,1.38e3,4.92e3), '8':(0.0959,0.0554,0.0125,0.196)}

#experimental rates from Berger et al., Cell Chemical Biology 2021,https://doi.org/10.1016/j.chembiol.2021.01.003  
#format: mean, SE in s^-1
experimental_rates={'1':(7.4e-2,2.5e-2),'2 (AMBER)':(75e-2,32e-2),'2 (CHARMM)':(75e-2,32e-2),'8':(25e-2,11e-2)}

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_xscale('log')
ax.set_yscale('log')

x=[]
xerr=[]
y=[]
yerr=[]
for ligand in computational_rates.keys():
	comp_rates=computational_rates[ligand]
	exp_rates=experimental_rates[ligand]
	#plot exp. rates on x-axis against comp. rates on y-axis
	x.append(exp_rates[0])
	xerr.append(exp_rates[1])
	y.append(comp_rates[0]) #all-data estimate
	yerr.append([comp_rates[0]-comp_rates[2],comp_rates[3]-comp_rates[0]]) #CI is not symmertic, need errors as difference
	ax.text(x=exp_rates[0],y=comp_rates[0],s=ligand,horizontalalignment='left',verticalalignment='bottom',fontsize=11)

	print(portray(comp_rates[0],comp_rates[2],comp_rates[3]))
yerr=np.array(yerr).transpose()
print(yerr)
ax.errorbar(x=x,xerr=xerr,y=y,yerr=yerr,c='k',linestyle='none')
#ax.set_xlim(left=1e-2,right=1e2)
#ax.set_ylim(bottom=1e-2,top=1e5)
#ticks=[10**x for x in range(-2,6)]
#ax.set_xticks(ticks)
#ax.set_yticks(ticks)
#ax.set_aspect('equal')

#lims = np.array([
    #np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
#    1/np.max([ax.get_xlim(), ax.get_ylim()]),
#    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
#])
lims=np.array([1e-5,1e5])
# now plot both limits against eachother
ax.plot(lims, lims, 'k-', alpha=1.0, zorder=0)
#represent 2 orders-of-magnitude goal
ax.fill_between(x=lims, y1=lims/1e2, y2=lims*1e2, alpha=0.2,color='k')
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
ticks=10**np.linspace(start=-4,stop=4,num=5)
print(ticks)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.grid()
ax.set_xlabel('experimental rate (s${}^{-1}$)')
ax.set_ylabel('computational rate (s${}^{-1}$)')
for fmt in ['pdf','png']:
	outfname='pyk2-allcpds-rate.{0}'.format(fmt)
	fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)

