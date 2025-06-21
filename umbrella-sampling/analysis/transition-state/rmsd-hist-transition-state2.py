#!/usr/bin/env python

import sys
import math
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pdb
scratch='/lustre/scratch/jmsc87/pyk2/umbrella-amber/'
dt=0.001 #in ns
stride=int(sys.argv[1])
binsize=float(sys.argv[2])

names=['pyk2-cpd1-amber','pyk2-cpd2-amber','pyk2-cpd8-amber']
labels={'pyk2-cpd1-amber':'ligand 1','pyk2-cpd2-amber':'ligand 2','pyk2-cpd8-amber':'ligand 8'}
states=['bound','transition']
linestyles={'bound':'-','transition':'--'}
colors={'pyk2-cpd1-amber':'C0','pyk2-cpd2-amber':'C1','pyk2-cpd8-amber':'C2'}
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
for iname,name in enumerate(names):
	#load everything
	psffname='prmtop/{0}.prmtop'.format(name)
	top=md.load_prmtop(psffname)
	refpdb='pdb/{0}-ref.pdb'.format(name)
	reftraj=md.load(refpdb)
	#residue numbers are zero based!
	nterm_backbone_atoms=top.select('(residue 0 to 84) and (name N or name CA or name C)')
	cterm_backbone_atoms=top.select('(residue 87 to 145 or residue 167 to 270) and (name N or name CA or name C)')
	for state in states:
		dcdfname='{0}/transition-state/{1}-{2}-state.dcd'.format(scratch,name,state)
		print('reading trajectory ',dcdfname,flush=True)
		traj=md.load(dcdfname,top=top,stride=stride)
		traj.superpose(reftraj,atom_indices=cterm_backbone_atoms)
		rmsd=md.rmsd(traj,reftraj,atom_indices=nterm_backbone_atoms)*10
		nbins=int(np.ceil(np.max(rmsd)/binsize))
		label=labels[name]+', '+state
		n,bins,patches=ax.hist(rmsd,bins=nbins,range=(0,nbins*binsize),histtype='step',density=True,
			color=colors[name],linestyle=linestyles[state],label=label)


ax.set_xlabel('N-terminal domain backbone RMSD (A)')
ax.set_ylabel('probability density')
ax.legend(loc='upper right',fontsize=8)
for fmt in ['png','pdf']:
	outfname='pyk2-allcpds-rmsd-hist.{0}'.format(fmt)
	fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=300)
plt.close(fig)

