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
#milestonefile='{0}-milestones-2'.format(name)
milestonefile=sys.argv[2]
initial_milestones=[]
input=open(milestonefile,'r')
for line in input:
	words=line.split()
	ianchor=int(words[1])
	janchor=int(words[2])
	initial_milestones.append((ianchor,janchor))
input.close()
#print(initial_milestones)

logfile='{0}-log'.format(name)
new_milestones=[]
input=open(logfile,'r')
for line in input:
	words=line.split()
	#these are the target anchors to which the 
	inewanchor=int(words[4])
	jnewanchor=int(words[5])
	t=(inewanchor,jnewanchor)
	if ((t not in initial_milestones) and (t not in new_milestones)):
		new_milestones.append(t)
input.close()
print(len(new_milestones),' new milestones discovered')
new_milestones=sorted(new_milestones)
print(new_milestones)


#put out a file with the new milestones
nmilestone=len(initial_milestones)
outfname='{0}-new-milestones'.format(name)
output=open(outfname,'w')
imilestone=nmilestone+1
for m in new_milestones:
	line='{0:d} {1:d} {2:d}\n'.format(imilestone,m[0],m[1])
	output.write(line)
	imilestone+=1
output.close()	

initial_milestones=np.array(initial_milestones)
new_milestones=np.array(new_milestones)


fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.scatter(initial_milestones[:,0],initial_milestones[:,1],color='black',label='initial milestones')
ax.scatter(new_milestones[:,0],new_milestones[:,1],color='red',label='new milestones')
ax.set_aspect('equal')
ax.set_xlabel('anchor number')
ax.set_ylabel('anchor number')
ax.legend(loc='lower right',fontsize=11)
outfname='{0}-milestones.png'.format(name)
fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
plt.close(fig)


