from scipy.stats import zscore 
import math
from random import randrange
import datetime,pandas
import numpy as np
from matplotlib import *
from math import *
import pylab,sys
from scipy.spatial.distance import mahalanobis as maha
import time 

def sqr(z):
	"""
	@param x:
	@return:
	"""
	return z * z

#m=int(sys.argv[1])		## Number of points of candidate
#q=int(sys.argv[2])		## Number of points of query
missing_file = sys.argv[1] 		## Number of points to be made None
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

def MahalanobisDist(x, y):
    covariance_xy = np.cov(x,y, rowvar=0)
    inv_covariance_xy = np.linalg.pinv(covariance_xy)

    xy_mean = np.mean(x),np.mean(y)
    x_diff = np.array([x_i - xy_mean[0] for x_i in x])
    y_diff = np.array([y_i - xy_mean[1] for y_i in y])
    diff_xy = np.transpose([x_diff, y_diff])
    
    md = []
    for i in range(len(diff_xy)):
        md.append(np.sqrt(np.dot(np.dot(np.transpose(diff_xy[i]),inv_covariance_xy),diff_xy[i])))
    return sum(md)

#------------------------------  NEW METHOD  ---------------------------------#
wincount=0
count=0
loc1=0
mahaldist1=100000000

start = time.time()
for i in range(len(finalval)-len(finalquery)+1):
	x=finalval[i:(i+len(finalquery))]
	prev=0
	j=0
	md=0
	count=0

	for k in x:
		if(k != None):
			count+=1
	nones=q-count
	
        if((float(nones)/q)>tol):
		print "SKIPPED: ",float(nones)/q
		continue
	wincount=0
	y=[]
	query=[]
	while(j<len(x)):
		if(j % winsize==0):
			prev=j
		if(wincount<winsize and x[j]!=None):
			wincount+=1
			y.append(x[j])
			query.append(finalquery[j])
			j+=1

		elif(x[j]==None):
			j+=1	
		elif(wincount==winsize):
			md+=MahalanobisDist(y,query)
			wincount=0
			j=prev
	
	md+=MahalanobisDist(y,query)
	ml=md/count
	if(mahaldist1>ml):
                mahaldist1=ml
                loc1=i


#-----------remove_end--------------#

end=time.time()
	
print("Location new: ",loc1)
print("Distance new: ",mahaldist1)

print("Time taken : ",end-start)
pylab.figure(1)
pylab.plot([x for x in range((loc1),(loc1+len(finalquery)))],finalquery,"-g",label="new query")
pylab.plot([x for x in range((loc1),(loc1+len(finalquery)))],finalval[(loc1):(loc1+len(finalquery))],"-b",label="new candidate")

pylab.legend(loc="upper right")
pylab.savefig("WindowVectors.png")



