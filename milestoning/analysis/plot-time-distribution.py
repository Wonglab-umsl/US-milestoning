#!/usr/bin/env python
import sys
import os
import numpy as np
import math
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt

name=sys.argv[1]
logfile='{0}-log'.format(name)
times=[]
input=open(logfile,'r')
for line in input:
	words=line.split()
	time=float(words[7])
	times.append(time)
input.close()
times=np.sort(np.array(times))
#print(times)
n=times.shape[0]
print(n)
y=np.linspace(0,1,num=n)
#print(y)


fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.plot(times,y)
ax.set_xlabel('stopping time (ns)')
ax.set_xscale('log')
ax.set_ylabel('cumulative distribution')
#ax.legend(loc='upper left',fontsize=11)
outfname='{0}-time-distribution.png'.format(name)
fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
plt.close(fig)

