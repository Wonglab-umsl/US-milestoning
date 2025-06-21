#!/usr/bin/env python
#must use python 2.7.16
import sys
import numpy as np
import math

fname=sys.argv[1]
#fname='rmsd-{0}'.format(name)
input=open(fname,'r')
#script will declare windows whose centers are at least interval A apart
interval=float(sys.argv[2]) 
#and are within limit A of the starting COM 
lowlimit=float(sys.argv[3])
highlimit=float(sys.argv[4])
#mcom=np.zeros(3)
#prev_com=None
window_centers=[]
new_window=False
lno=1
for line in input:
	words=line.split()
	#time, backbone RMSD, com x, com y, com z
	iwindow=int(words[0])
	#time=float(words[2])
	comx=float(words[1])
	comy=float(words[2])
	comz=float(words[3])
	com=np.array([comx,comy,comz])
	if (lno==1):
		com0=com
	#the new window has to be at least interval A away from ANY previous window
	new_window=True
	mind=1.0e20
	for w in window_centers:
		disp=com-w
		d=math.sqrt(np.dot(disp,disp))
		if (d<mind):
			mind=d
	disp2=com-com0
	d2=math.sqrt(np.dot(disp2,disp2))
	if ((mind>interval) and (d2>lowlimit) and (d2<highlimit)):
		window_centers.append(com)
		print iwindow,com[0],com[1],com[2],mind,d2
	lno=lno+1	
input.close()

