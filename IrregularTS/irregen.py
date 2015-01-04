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


#m=int(sys.argv[1])	## Number of points of candidate
#q=int(sys.argv[2])	## Number of points of query
n=int(sys.argv[1])	## Number of points to be made None

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
'''
sampledata=val[756562-256:756562+256]
#print(sampledata)
f=open("SampleData.txt","w")
for i in sampledata:
	f.write(str(i)+" ")
f.close()
f=open("SampleQuery.txt","w")
samplequery=val[756562:756562+129]
for i in samplequery:
        f.write(str(i)+" ")
f.close()

print("DOne writing")
'''
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
'''
df= pandas.DataFrame({"finalquery":query})
df['finalquery'] = (df.finalquery - df.finalquery.mean())/df.finalquery.std(ddof=0)

df= pandas.DataFrame({"val":val})
df['val'] = (df.val - df.val.mean())/df.val.std(ddof=0)

val=df['val']
'''
ex = sum(query)
ex1 = sum(map(sqr, query))
mean = ex / len(query)
std = ex1 / len(query)
std = math.sqrt(std - mean * mean)
finalquery = [(X - mean) / std for X in query]


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
		print("None loc : ",i)
#print finalval

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

count=0
loc=0
mahaldist=100000000
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
	if(mahaldist>ml):
                mahaldist=ml
                loc=i
	#print(mahaldist,ml,i)
	
print("Location:",loc)
print("Distance:",mahaldist)
pylab.plot([x for x in range((loc),(loc+len(finalquery)))],finalval[(loc):(loc+len(finalquery))],"-b",label="candidate")	
pylab.plot([x for x in range((loc),(loc+len(finalquery)))],finalquery,"-g",label="query")
pylab.legend(loc="upper right")
pylab.savefig("irregplot.png")

#with open("results2","a") as f:
#	f.write("{0},{1},{2}\n".format(n,loc,mahaldist))
#md=MahalanobisDist(finalval,finalquery)
#print(md)

#print finalval
