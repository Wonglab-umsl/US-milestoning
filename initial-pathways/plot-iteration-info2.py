#!/usr/bin/env python3.6

#must use python 2.7.16
import sys
import os
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 18,'legend.fontsize':12})
import matplotlib.pyplot as plt

#cpds=['cpd2','cpd3','cpd5','cpd6','cpd8','cpd11','cpd12','cpd12n']
cpds=['cpd2','cpd8']

fig=plt.figure()
ax=fig.add_subplot(1,1,1)

for cpd in cpds:
	name='pyk2-{0}-amber'.format(cpd)
	print(name)
	infname='data/{0}-iteration-info'.format(name)
	input=open(infname,'r')
	iterlist=[]
	lowlist=[]
	highlist=[]
	for line in input:
		#if (line.find("ligand COM")>=0
		words=line.split()
		iterlist.append(int(words[0]))
		lowlist.append(float(words[1]))
		highlist.append(float(words[2]))

	#ax.plot(iterlist,lowlist)
	label=cpd.replace("cpd","")
	ax.fill_between(iterlist,lowlist,highlist,alpha=0.5,label=label)

#ax.plot(iterlist,highlist)
ax.set_ylim(bottom=0,top=20)
ax.legend(loc='upper right',ncol=2)
for fmt in ['png','pdf']:
	outfname='pyk2-allcpds-amber-iteration-info.{0}'.format(fmt)
	#if (fmt=='pdf'):
	ax.set_xlabel('iteration number')
	ax.set_ylabel('ligand COM distance (A)')
	fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)
plt.close(fig)
