import sys
import re
import lipid as lib
import numpy as np
from math import sqrt
from math import pi

#"preprocessor"
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    """Small modification such that first and last chunk are of same size. Middle chunk is messed up"""
    size = len(l)//n
    lst = [size for i in range(n)]
    leftover = len(l) - size * n
    mid = n//2
    
    if leftover < size//2:     
        lst[mid] += leftover
    else:
        lst.insert(mid+1,leftover)
    """
    tot = 0
    for element in lst:
        tot += element
        
    if tot != len(l):
        return "Failure"
    """    
    curr = 0
    for i in range(len(lst)):
        yield l[curr:curr+lst[i]]
        curr += lst[i]
		
NDIM = 3
DR = 0.02
dt = 0.02

#Arguments of script

lipidFileName = "dr_lipid_top-0-45.dat"
cholFileName = "dr_chol_top-0-45.dat"

lipidarr = []
cholarr = []
r = np.empty([NDIM-1])

file = open(lipidFileName,'r')
for line in file:
    x,y,s = line.split()
    lipidarr.append([x,y,s])
file.close()

file = open(cholFileName)
for line in file:
    x,y,s = line.split()
    cholarr.append([x,y,s])
file.close()
	
lipidarr = sorted(lipidarr, key = lambda x : x[2])
cholarr = sorted(cholarr, key = lambda x : x[2])

Nlipid = len(lipidarr)
Nchol = len(cholarr)

lipidChunk = chunks(lipidarr,Nlipid//20)
del lipidarr

chunkCount = 0

for chunk in lipidChunk:
    for i in range(len(chunk)):
        for j in range(1,len(chunk)):
            for k in range(NDIM):
                r[k] = chunk[i][k] - chunk[j][k]
                r[k] = lib.periodic(r[k],L[k])
            r2 = sqrt(r[0]*r[0] + r[1]*r[1])
        
    chunkCount += 1


    




