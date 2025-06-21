#!/usr/bin/env python
import sys
import numpy as np
import math
import os

name=sys.argv[1]
#nunbiased=int(sys.argv[2])
nwindow=int(sys.argv[2])
#nsim=nunbiased+nwindow
nanchor=int(sys.argv[3])
anchorsfname=sys.argv[4]
outfname=sys.argv[5]
dir=os.environ['HOME']+'/pyk2/umbrella-amber/anal3/data-rel-to-first'
scratch='/mnt/stor/ceph/scratch/jmsc87/pyk2/'
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
nanchor=anchors.shape[0]
#print(nanchor)
dt=0.001
#com_data={}
window_and_time=[]
com_all=[]
for iwindow in range(0,nwindow):
	fname='{0}/{1}-tramd-com-{2:d}'.format(dir,name,iwindow+1)
	input=open(fname,'r')
	times=[]
	com=[]
	for line in input:
		words=line.split()
		#frame indices for catdcd are zero based
		times.append(int(float(words[0])/dt)-1)
		com.append(np.array([float(words[1]),float(words[2]),float(words[3])]))
	input.close()
	times=np.array(times)
	tmp=np.full_like(times,iwindow+1,dtype=np.int32)
	times=np.vstack([tmp,times]).transpose()
	window_and_time.append(times)
	#print(times)
	com=np.array(com)
	print(iwindow,com.shape)
	com_all.append(com)

times=np.vstack(window_and_time)
print(times.shape)
print(times)
com=np.vstack(com_all)
print(com.shape)

nframes=com.shape[0]
dist=np.full((nframes,nanchor),np.nan,dtype=np.float64)
for ianchor in range(0,nanchor):
	print("computing distances from anchor ",ianchor+1,flush=True)
	diff=com[:,:]-np.expand_dims(anchors[ianchor,:],axis=0)
	dist[:,ianchor]=np.sqrt(np.einsum('ij,ij->i',diff,diff)) #intended shape (nframes)

#choose the two anchors that are closest for each frame
info=np.argsort(dist,axis=1)
milestones=np.sort(info[:,0:2],axis=1) #so that the lower-numbered milestone is first
print(milestones.shape)
print(milestones)
#sys.exit(-1)

#fname='{0}/umbrella-amber3/tram-anal2b/{1}-ligand-com-data.npy'.format(scratch,name)
#print('loading ligand COM data at ',datetime.now(),flush=True)
#this has been converted to float32 -- need to convert it back (better would be to write it as float64,
#but not sure how to control the type of get_output above.
#ligand_com=np.load(fname).astype(np.float64)
#ligand_com=load_numpy_array_showing_progress(fname,1,'loading ligand COM data')
#ligand_com=np.load(fname,allow_pickle=True)
#ligand_com=[ligand_com[i].astype(np.float64) for i in range(0,len(ligand_com))]
#ligand_com=ligand_com.astype(np.float64)
#ligand_com=list(ligand_com)
#print(len(ligand_com),ligand_com[0].shape,ligand_com[0].dtype)
        #print(ligand_com.shape)
#ligand_com_unbiased=ligand_com[0:nunbiased]
#ligand_com_umbrella=ligand_com[nunbiased:nsim]

mcenter=np.full(3,np.nan,dtype=np.float64)
imilestone=1
output=open(outfname,'w')
for i in range(0,nanchor-1):
	ianchor=i
	janchor=i+1
	#skip over anchor point 47 for ligand 8 to try to get something that is easier to select starting frames for
	if ((name=='pyk2-cpd8-amber') and (ianchor+1==46)):
		janchor=i+2
	if ((name=='pyk2-cpd8-amber') and (ianchor+1==47)):
		continue
	if ((name=='pyk2-cpd10-amber') and (ianchor+1==50)): #milestone 50,51
		janchor=i+2
	if ((name=='pyk2-cpd10-amber') and (ianchor+1==51)):
		continue
	#select those frames which are on milestone
	#this_anchor=np.array([ianchor,janchor])
	mask=np.all([(milestones[:,0]==ianchor),(milestones[:,1]==janchor)],axis=0)
	indices=np.where(mask)[0] #get indices in which mask is true)
	#print(mask)
	#print(indices)
	#print(milestones[mask,:])
	count=np.count_nonzero(mask)
	print("milestone {0:d} {1:d} has {2:d} frames".format(ianchor+1,janchor+1,count),flush=True)
	if (count==0):
		#can't find a frame
		#line='{0:d} {1:d} {2:d}\n'.format(imilestone,ianchor+1,janchor+1)
		#output.write(line)
		print("warning: couldn't find any frames on milestone, relaxing criteria",flush=True)
		#try allowing frames where only one of the closest 
		mask=np.any([(milestones[:,0]==ianchor),(milestones[:,1]==janchor)],axis=0)
		indices=np.where(mask)[0]
		count=np.count_nonzero(mask)
		print("now have {2:d} frames".format(ianchor+1,janchor+1,count),flush=True)
		#continue
	#sys.exit(-1)
	#continue
	#milestone between i and i + 1
	mcenter[:]=0.5*(anchors[ianchor,:]+anchors[janchor,:])
	#determine distance to the milestone midpoint, choosing only from those frames that are on milestone
	diff=com[mask,:]-np.expand_dims(mcenter[:],axis=0) #shape (nframes,3) 
	dist=np.sqrt(np.einsum('ij,ij->i',diff,diff)) #intended shape (nframes)
	#print(dist.shape)
	#sys.exit(-1)
	#the umbrella with minimum distance
	ii=np.argmin(dist) #this is out of the frames in mask
	mindist=dist[ii] 
	idx=indices[ii] #this is an index in the entire frame set
	print(milestones[idx,:])
	line='{0:d} {1:d} {2:d} {3:d} {4:d} {5:.6f} {6:.6f} {7:.6f} {8:.6f}\n'.format(
		imilestone,ianchor+1,janchor+1,times[idx,0],times[idx,1],com[idx,0],com[idx,1],com[idx,2],mindist)
	output.write(line)
	imilestone+=1
output.close()
