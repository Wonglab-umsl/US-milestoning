#!/usr/bin/env python

#must use python 2.7.16
import sys
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 18,'legend.fontsize':12})
import matplotlib.pyplot as plt
from scipy import stats


barrier_heights={'1':12.75,'10':17.5,'2':7.75,'6':9.25,'8':15.25}

#exp_rates={'1':7.4e-2,'2':75e-2,'10':4.0e-2,'6':19e-2,'8':25e-2}
#first number is mean, second is sd from table 1 in  experimental rate paper
exp_rates={'1':(7.4e-2,2.5e-2),'2':(75e-2,32e-2),'10':(4.0e-2,1.2e-2),'6':(19e-2,7.8e-2),'8':(25e-2,11e-2)}

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
x=[]
y=[]
xerr=[]
for cpd in barrier_heights.keys():
	x.append(exp_rates[cpd][0])
	xerr.append(exp_rates[cpd][1])
	y.append(barrier_heights[cpd])
#print(x)
#print(y)
x=np.array(x)
y=np.array(y)
#ax.scatter(x,y)
ax.errorbar(x,y,xerr=xerr,marker='+',ls='none',capsize=6)
for i,cpd in enumerate(barrier_heights.keys()):
	ax.annotate(cpd,(x[i],y[i]))

slope, intercept, r, p, stderr = stats.linregress(np.log(x),y)
print(slope,intercept,r,p,stderr)
ybar=slope*np.log(x)+intercept
ax.plot(x,ybar)

ax.set_xlabel('experimental dissociation rate (s${}^{-1}$)')
ax.set_ylabel('barrier height (kcal/mol)')
ax.set_xscale('log')
ax.set_xlim(left=None,right=1e-0)
for fmt in ['png','pdf']:
	output='pyk2-allcpds-barrier-height-vs-rate.{0}'.format(fmt)
	fig.savefig(output,format=fmt,bbox_inches='tight',dpi=600)

plt.close(fig)
