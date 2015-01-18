from random import randrange
import numpy as np
from matplotlib import *
import pylab,sys

missing_file=sys.argv[1] 		## Number of points to be made None
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

	pmm = abs(P-mean)
	md=(float(pmm)/float(std))
	return md

#------------------------------  NEW METHOD  ---------------------------------#
wincount=0
count=0
loc1=0
mahaldist1=10000000000000

import time
starttime = time.time()

for i in range(len(finalval)-len(finalquery)+1):
	if( i in range(5)):
		x=finalval[i:i+len(finalquery)]
	else:
		x=finalval[i-5:(i+len(finalquery))]		# chunk of candidate of size query
	prev=0
	j=0
	count=0
	for k in range(len(x)):
		if(x[k] != None):
			count+=1

	nones=q-count
	if((float(nones)/q)>tol):
		print "SKIPPED: ",float(nones)/q
		continue
	md=0
	md_temp=0
	for j in range(5):
		md_temp=incr_maha(x[0:10],finalquery[j])
		if(not(md_temp == 10000000000000)):
			md+=md_temp
	start = 5
	end = 15
	for j in range(5,len(finalquery)-4):
		md_temp=incr_maha(x[start:end],finalquery[j])
		start+=1
		end+=1
		if(not(md_temp == 10000000000000)):
			md+=md_temp
	for j in range(len(finalquery)-4,len(finalquery)):
		md_temp=incr_maha(x[start:end],finalquery[j])
		if(not(md_temp == 10000000000000)):
			md+=md_temp
	if(md < mahaldist1):
		mahaldist1=md
		loc1=i

endtime = time.time()
print("Location new incr: ",loc1)
print("Distance new incr: ",mahaldist1)
print("TIme taken : ",endtime-starttime)

pylab.figure(1)
pylab.plot([x for x in range((loc1),(loc1+len(finalquery)))],finalquery,"-g",label="new query")
pylab.plot([x for x in range((loc1),(loc1+len(finalquery)))],finalval[(loc1):(loc1+len(finalquery))],"-b",label="new candidate")

pylab.legend(loc="upper right")
pylab.savefig("P2D.png")
