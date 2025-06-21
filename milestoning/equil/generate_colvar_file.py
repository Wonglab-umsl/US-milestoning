#!/usr/bin/env python
import sys
import numpy as np
import math


#        distance {
#                group1 {
#                        atomsFile :refpdb;
#                        atomsCol O
#                        centerReference yes
#                        rotateReference yes
#                        refPositionsFile :refpdb;
#                        refPositionsCol B
#                        refPositionsGroup {
#                                atomsFile :refpdb;
#                                atomsCol B
#                        }
#                }
#                group2 {
#                        dummyAtom (:com;)
#                }
#        }

#write a distance group.   This defines the distance between the ligand COM and the specified point
#after RMS alignment with the protein backbone

def write_colvar_distance(output,coeff,refpdb,point):
	output.write('\tdistance {\n')
	output.write('\t\tcomponentCoeff {0:.1f}\n'.format(coeff))
	output.write('\t\tgroup1 {\n')
	output.write('\t\t\tatomsFile {0}\n'.format(refpdb))
	output.write('\t\t\tatomsCol O\n')
	output.write('\t\t\tcenterReference yes\n')
	output.write('\t\t\trotateReference yes\n')
	output.write('\t\t\trefPositionsFile {0}\n'.format(refpdb))
	output.write('\t\t\trefPositionsCol B\n')
	output.write('\t\t\tfittingGroup {\n')
	output.write('\t\t\t\tatomsFile {0}\n'.format(refpdb))
	output.write('\t\t\t\tatomsCol B\n')
	output.write('\t\t\t}\n')
	output.write('\t\t}\n')
	output.write('\t\tgroup2 {\n')
	output.write('\t\t\tdummyAtom ({0:.6f},{1:.6f},{2:.6f})\n'.format(point[0],point[1],point[2]))
	output.write('\t\t}\n')
	output.write('\t}\n')




nanchor=int(sys.argv[1])
anchorsfname=sys.argv[2]
milestonesfname=sys.argv[3]
this_milestone_index=int(sys.argv[4])
refpdb=sys.argv[5]
kharm=float(sys.argv[6])
kharm_wall=float(sys.argv[7])
outfname=sys.argv[8]


#format: ianchor x y z energy, assumed in order
input=open(anchorsfname,'r')
anchors=np.full((nanchor,3),np.nan,dtype=np.float64)
for line in input:
	words=line.split()
	ianchor=int(words[0])-1
	for k in range(0,3):
		anchors[ianchor,k]=float(words[k+1])
input.close()
anchors=np.array(anchors)
#print(anchors)
#format: imilestone ianchor janchor startng-umbrella [the last is not needed]
milestones={}
input=open(milestonesfname,'r')
this_milestone=None
for line in input:
	words=line.split()
	#these indices are 1-based
	imilestone=int(words[0])
	anchor=(int(words[1]),int(words[2]))
	if (imilestone==this_milestone_index):
		this_milestone=anchor
	if (imilestone in milestones):
		print("error: duplicate milestone")
		sys.exit(-1)
	milestones[imilestone]=anchor
#print(milestones)
input.close()

#sys.exit(-1)
#print(this_milestone)
if (this_milestone is None):
	print("error: can't find specified milestone")
	sys.exit(-1)

output=open(outfname,'w')
output.write('# generated using anchors file {0} and milestones file {1}\n'.format(
	anchorsfname,milestonesfname))
output.write('# for milestone number {0:d} anchors {1}\n'.format(
	this_milestone_index,this_milestone))
output.write('\n')

#write this milestone's colvar definition
ianchor=this_milestone[0]
ianchor_pos=anchors[ianchor-1,:]
janchor=this_milestone[1]
janchor_pos=anchors[janchor-1,:]

#print(ianchor,janchor)
output.write('colvar {\n')
output.write('\tname neighbor\n')
write_colvar_distance(output,1.0,refpdb,ianchor_pos)
write_colvar_distance(output,-1.0,refpdb,janchor_pos)
output.write('}\n')
output.write('\n')

#write the harmonic potential on this distance difference
output.write('harmonic {\n')
output.write('\tcolvars neighbor\n')
output.write('\tcenters 0\n')
output.write('\tforceConstant {0:.1f}\n'.format(kharm))
output.write('}\n')
output.write('\n')


#sys.exit(-1)
#now write all the repulsive potentials re other milestones
for milestone in milestones.values():
	#we need only those milestones which share at least one anchor with this_milestone
	#each milestone should be a tuple of the form (i,j) with i<j
	#convert them to sets for comparison, just in case we have milestones (i,j) and (j,i)
	if ((set(milestone)==set(this_milestone)) or set(this_milestone).isdisjoint(set(milestone))):
		continue
	elif (milestone[0] in this_milestone):
		this_anchor=milestone[0]
		other_anchor=milestone[1]
	elif (milestone[1] in this_milestone):
		this_anchor=milestone[1]
		other_anchor=milestone[0]
	else:
		print("this shouldn't happen")
		sys.exit(-1)
	#this_anchor is the one that overlaps
	#other_anchor is the one that doesn't
	this_anchor_pos=anchors[this_anchor-1,:]
	other_anchor_pos=anchors[other_anchor-1,:]
	output.write('colvar {\n')
	#other_anchor gets a weight of +1, this_anchor 
	name='{0:d}_{1:d}'.format(other_anchor,this_anchor)
	output.write('\tname {0}\n'.format(name))
	write_colvar_distance(output,1.0,refpdb,other_anchor_pos)
	write_colvar_distance(output,-1.0,refpdb,this_anchor_pos)
	output.write('}\n')
	output.write('\n')
	#the new format for colvar files requires the harmonic part to be written separately
	output.write('harmonicWalls {\n')
	output.write('\tcolvars {0}\n'.format(name))
	#I don't know why Samith's scripts have 0.1 here
	output.write('\tlowerWalls 0\n')
	output.write('\tlowerWallConstant {0:f}\n'.format(kharm_wall))
	output.write('}\n')
	output.write('\n')

