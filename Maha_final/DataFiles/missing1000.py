from random import randrange
import sys

n=int(sys.argv[1])
setter=set()
while(len(setter)<n):
	r=randrange(1,6000)
	setter.add(r)

f=open("missing"+str(n),"w")
for i in setter:
	f.write(str(i)+"\n")

f.close()

