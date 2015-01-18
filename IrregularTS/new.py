from scipy.stats import zscore 
import math
from random import randrange
import datetime,pandas
import numpy as np
from matplotlib import *
from math import *
import pylab,sys

def sqr(z):
	"""
	@param x:
	@return:
	"""
	return z * z

#m=int(sys.argv[1])		## Number of points of candidate
#q=int(sys.argv[2])		## Number of points of query
n=int(sys.argv[1]) 		## Number of points to be made None
tol=float(sys.argv[2])		## Tolerance level


val=[]
query=[]

f=open("SampleData.txt","r")
s1=f.read()
l1=s1.split()
for x1 in l1:
	if(x1=="None"):
		val.append(None)
	else:
		y1=float(x1)
		val.append(y1)
f.close()
print("Number of points in Candidate : %d"%len(val))

print("Total points= ",len(val))
f=open("SampleQuery.txt","r")
s2=f.read()
l2=s2.split()
for x2 in l2:
	if(x2=="None"):
		query.append(None)
	else:
		y2=float(x2)
		query.append(y2)
f.close()

d1=dict()
#d2=dict()

ex = sum(query)
ex1 = sum(map(sqr, query))
mean = ex / len(query)
std = ex1 / len(query)
std = math.sqrt(std - mean * mean)
finalquery = [(X - mean) / std for X in query]

q=len(finalquery)
#winsize=int(((1-tol)/2)*q)

ex = sum(val)
ex1 = sum(map(sqr, val))
mean = ex / len(val)
std = ex1 / len(val)
std = math.sqrt(std - mean * mean)
val = [(X - mean) / std for X in val]


for i in range(len(val)):
	d1[i+1]=val[i]

for i in range(n):
	r=randrange(1,len(val))
	d1[r]=None


finalval=[]


for i in range(1,len(val)+1):
	finalval.append(d1[i])

for i in range(len(finalval)):
	if(finalval[i]==None):
		print("None loc : %d"%i)
#print finalval

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

#------------------------------  RUNNING MAHA  ---------------------------------#
wincount=0
count=0
loc1=0
mahaldist1=100000000


for i in range(len(finalval)-len(finalquery)+1):
	x=finalval[i:(i+len(finalquery))]
	prev=0
	j=0
	md=0
	count=0
#	print(x)
	for k in x:
		if(k != None):
			count+=1
	
	while(j<len(x)):
		if(x[j]==None):
			y=x[prev:(j)]
			query=finalquery[prev:(j)]
			j+=1
			prev=j
			if(len(y)>0):
				md+=MahalanobisDist(y,query)
		else:
		
			j+=1
	y=x[prev:j]
	query=finalquery[prev:j]
	if(len(y)>0):
		md+=MahalanobisDist(y,query)
	ml=md/count
	if(mahaldist1>ml):
                mahaldist1=ml
                loc1=i
	#print(mahaldist,ml,i)

	
print("Location regular: ",loc1)
print("Distance regular: ",mahaldist1)

pylab.plot([x for x in range((loc1),(loc1+len(finalquery)))],finalval[(loc1):(loc1+len(finalquery))],"-b",label="original candidate")
pylab.plot([x for x in range((loc1),(loc1+len(finalquery)))],finalquery,"-g",label="query")
pylab.legend(loc="upper right")
pylab.savefig("irregplot.png")

