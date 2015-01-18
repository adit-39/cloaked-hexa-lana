from random import randrange
import numpy as np
from matplotlib import *
import pylab,sys

missing_file=sys.argv[1] 		## Number of points to be made Non
tol=float(sys.argv[2])		## Tolerance level


val=[]
finalquery=[]

f=open("out","r")
s1=f.read()
l1=s1.split()
for x1 in l1:
	if(x1=="None"):
		val.append(None)
	else:
		y1=float(x1)
		val.append(y1)
f.close()

print("Total points= ",len(val))
f=open("query","r")
s2=f.read()
l2=s2.split()
for x2 in l2:
	if(x2=="None"):
		finalquery.append(None)
	else:
		y2=float(x2)
		finalquery.append(y2)
f.close()

d1=dict()



q=len(finalquery)
winsize=int(((1-tol)/2)*q)

for i in range(len(val)):
	d1[i+1]=val[i]

# Making points missing
f = open(missing_file,"r")
miss = f.read().split("\n")
for i in miss:
	d1[i]=None
f.close()

finalval=[]

for i in range(1,len(val)+1):
	finalval.append(d1[i])

#------------------------------  MAHALANOBIS  ----------------------------------#
def covariance(x, y):
    covariance_xy = np.cov(x,y, rowvar=0)
    inv_covariance_xy = np.linalg.pinv(covariance_xy)
    return inv_covariance_xy

def MahalanobisDist(x, y):
    inv_covariance_xy = covariance(x,y)
    xy_mean = np.mean(x),np.mean(y)
    x_diff = np.array([x_i - xy_mean[0] for x_i in x])
    y_diff = np.array([y_i - xy_mean[1] for y_i in y])
    diff_xy = np.transpose([x_diff, y_diff])
    
    md = []
    for i in range(len(diff_xy)):
        md.append(np.sqrt(np.dot(np.dot(np.transpose(diff_xy[i]),inv_covariance_xy),diff_xy[i])))
    return sum(md)

def incr_maha(D, P):
	D = [x for x in D if x is not None]
	if(len(D)==0):
		return 10000000000000
	mean = np.mean(D)
	std = float((sum([(x_i - mean)**2 for x_i in D])/len(D))**0.5)
#	print("std %f"%std)
	if(std==0):
		return 10000000000000
	pmm = abs(P-mean)
	md=float(pmm)/float(std)
	return md

#------------------------------  NEW METHOD  ---------------------------------#
wincount=0
count=0
loc1=0
mahaldist1=10000000000000

import time
start = time.time()

for i in range(len(finalval)-len(finalquery)+1):
	x=finalval[i:(i+len(finalquery))]		# chunk of candidate of size query
	prev=0
	j=0
	count=0
	c_y=[]
	c_query=[]
	for k in range(len(x)):
		if(x[k] != None):
			count+=1
			c_y.append(x[k])
			c_query.append(finalquery[k])

	nones=q-count
	
	if((float(nones)/q)>tol):
		print "SKIPPED: ",float(nones)/q
		continue
	
	inv_cov=covariance(c_y,c_query)
	wincount=0
	md=0
	
	while(True):
		md_temp=0
		if(wincount==2):
			win_val = x[prev:]
			win_q= finalquery[prev:]
		else:
			win_val = x[prev:prev+winsize]
			win_q = finalquery[prev:prev+winsize]
		xx=[]
		yy=[]
		for j in range(len(win_q)):		
			if(j==0 or j==1 or j==(len(win_q)-1)):
				if(win_val[j] is None):
					continue
				md+= ((win_val[j]-win_q[j])**2)**0.5
				xx.append(win_val[j])
				yy.append(win_q[j])
			else:
				xx.append(win_val[j])
				yy.append(win_q[j])				



		for j in range(1,len(xx)-1):
			mdimo = 10000000000000
			mdipo = 10000000000000
			mdi = 10000000000000
		
			templist=xx[:j-1]
			templist.append(xx[j-1])
			mdimo = incr_maha(templist,win_q[j])
				
			templist.pop()
			templist.append(xx[j])					
			mdi = incr_maha(templist,win_q[j])

			templist.pop()
			templist.append(xx[j+1])	
			mdipo = incr_maha(templist,win_q[j])

			md_temp = min(mdimo,mdi,mdipo)

		if(md_temp==10000000000000):
			pass
		else:
			md+=md_temp
		if(wincount==2):
			break
		wincount+= 1
		prev+=winsize
	if(md < mahaldist1):
		mahaldist1=md
		loc1=i

end=time.time()

print("Location new incr: ",loc1)
print("Distance new incr: ",mahaldist1)
print("Time taken : ",end-start)
pylab.figure(1)
pylab.plot([x for x in range((loc1),(loc1+len(finalquery)))],finalquery,"-g",label="new query")
pylab.plot([x for x in range((loc1),(loc1+len(finalquery)))],finalval[(loc1):(loc1+len(finalquery))],"-b",label="new candidate")

pylab.legend(loc="upper right")
pylab.savefig("P2DWarped.png")
