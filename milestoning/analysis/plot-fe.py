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
#pathwayfile='{0}-pathway-4'.format(name)
pathwayfile=sys.argv[2]
input=open(pathwayfile,'r')
anchors=[]
fe_from_umbrella=[]
for line in input:
	words=line.split()
	anchors.append(int(words[0]))
	fe_from_umbrella.append(float(words[4]))
input.close()
print(anchors)
fe_from_umbrella=np.array(fe_from_umbrella)
print(fe_from_umbrella)
renormalize=True
if (renormalize):
	kt=1.987e-3*300
	p=np.exp(-fe_from_umbrella/kt)
	p[:]=p[:]/np.sum(p)
	fe_from_umbrella=-kt*np.log(p)
print(fe_from_umbrella)
	
otherfile='{0}-fe-from-milestoning'.format(name)
input=open(otherfile,'r')
halfway=[]
fe_from_milestoning=[]
fe_ci_lo=[]
fe_ci_hi=[]
#std_fe_from_milestoning=[]
for line in input:
	words=line.split()
	#halfway.append(float(words[2]))
	ianchor=int(words[1])
	janchor=int(words[2])
	#use only "neighboring" anchors on the initial pathway
	if ((janchor-ianchor)==1):
		halfway.append((ianchor+janchor)/2)
		#as currently written, this is the free energy from the "all data" estimate
		fe_from_milestoning.append(float(words[3]))
		#these are the 68% CI
		fe_ci_lo.append(float(words[4]))
		fe_ci_hi.append(float(words[5]))
	#if (len(words)>4):
		#std_fe_from_milestoning.append(float(words[4]))
input.close()
print(halfway)
fe_from_milestoning=np.array(fe_from_milestoning)
fe_ci_lo=np.array(fe_ci_lo)
fe_ci_hi=np.array(fe_ci_hi)
#std_fe_from_milestoning=np.array(std_fe_from_milestoning)
print(fe_from_milestoning)
print(fe_ci_lo)
print(fe_ci_hi)
width=0.5*(fe_ci_hi[:]-fe_ci_lo[:])
rms_width=np.sqrt(np.sum(width[:]*width[:])/float(len(width)))
print("RMS half width of free energy conf interval: {0:.2f} kcal/mol".format(rms_width))
print("max half width of free energy conf interval: {0:.2f} kcal/mol".format(np.max(width)))

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.plot(anchors,fe_from_umbrella,label='from umbrella sampling')
ax.plot(halfway,fe_from_milestoning,label='from milestoning')
#if (len(std_fe_from_milestoning)>0):
ax.fill_between(halfway,fe_ci_lo,fe_ci_hi,alpha=0.2,color='C1')

ax.legend(loc='lower right',fontsize=11)
ax.set_xlabel('anchor number')
ax.set_ylabel('free energy (kcal/mol)')
for fmt in ['png','pdf']:
	outfname='{0}-fe.{1}'.format(name,fmt)
	fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)
plt.close(fig)
