import re
import sys

import numpy as np
cimport numpy as np

from math import sqrt
from bisect import bisect_left

def cholConc(filename):
    file = open(filename,"r")
    for line in file:
        m = re.match("CHOL", line)
        if m:
            line = line.split()
            file.close()
            return int(line[1]) * 2

cpdef processTraj(str trajFileName,int Nchol,int NDIM,int Nconf):
    trajFile = open(trajFileName,'r')
    cdef int N = int(trajFile.readline().split()[0])
    cdef int Nperlipid = 12
    cdef int Nperchol = 8
    cdef int Nlipids = (N - Nperchol * Nchol) // Nperlipid
    cdef int Ncholbeads = Nchol * Nperchol
    cdef int Nlipidbeads = Nlipids * Nperlipid
    
    #reading traj file
    trajFile.close()
    trajFile = open(trajFileName,'r')

   
    cdef np.ndarray[np.double_t, ndim=2] L
    cdef np.ndarray[np.double_t, ndim=3] x,y
    
    L = np.zeros((NDIM, Nconf), dtype=np.double)
    x = np.zeros((Nlipidbeads,NDIM,Nconf), dtype=np.double)
    y = np.zeros((Ncholbeads,NDIM,Nconf), dtype=np.double)
    
    cdef int t,i,j,k,index
    #cdef
       
    for t in range(Nconf):
        trajFile.readline()
        
        #Box Sizes
        for k in range(NDIM):
            L[k,t] = float(trajFile.readline().strip())
        
        #Bead Coordinates
        for i in range(Nlipidbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
                
            for k in range(NDIM):
                x[i,k,t] = float(config[k])
            
        for i in range(Ncholbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
            
            for k in range(NDIM):
                y[i,k,t] = float(config[k])
            
        for i in range(Nlipidbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
             
            for k in range(NDIM):
                index = i+Nlipidbeads//2
                x[index,k,t] = float(config[k])
        
        for i in range(Ncholbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
            
            for k in range(NDIM):
                index = i+Ncholbeads//2
                y[index,k,t] = float(config[k])

    trajFile.close()

    return N,L,x,y
    
cpdef comassing(np.ndarray[np.double_t, ndim=3] ar,str bead, np.ndarray[np.double_t, ndim=2] L):
    cpdef int size 
    if bead == "DPPC": 
        size= 12
        #center = 2
        #3rd bead but zero indexed
    elif bead == "CHOL":
        size = 8
        #center = 4
    else:
        print("bead type error")

    cdef int N = len(ar)//size
    cdef int NDIM = len(ar[0])
    cdef int Nconf = len(ar[0][0])
    
    cdef np.ndarray[np.double_t, ndim=3] comAr
    comAr = np.zeros((N,NDIM,Nconf), dtype=np.double)
    
    cdef int t,i,ii,k,j
    
    for t in range(Nconf):
        for k in range(NDIM):
            for i in range(N):
                ii = i*size
                
                #comAr[i,k,t] = ar[ii+center,k,t]
                
                for j in range(size): 
                    ar[ii+j,k,t] = periodic_p(ar[ii+j,k,t],ar[ii,k,t],L[k,t])
                    #fix your boundary condition
                    #look at structure factor
                    comAr[i,k,t] += ar[ii+j,k,t]/size
  
    return comAr

def periodic_p(current,first,L): #for periodic with reference to head bead
    if current - first > 0.5*L:
        return current - L
    elif current - first < -0.5*L:
        return current + L
    else:
        return current

def periodic(dr,L): #give the displacement and box size. Both flaots
    if dr > 0.5*L:
        return dr - L
    elif dr < -0.5*L:
        return dr + L
    else:
        return dr

cpdef translateZ(np.ndarray[np.double_t, ndim=3] x,np.ndarray[np.double_t, ndim=3] y):
    cdef int Nconf,Nlipidbeads,Ncholbeads,N
    cdef double z_avg
    cdef int t,i
    
    Nconf = len(x[0][0])
    Nlipidbeads = len(x)
    Ncholbeads = len(y)
    N = Nlipidbeads + Ncholbeads

    for t in range(Nconf):
        z_avg = 0.
            
        for i in range(Nlipidbeads):
            z_avg += x[i,2,t]

        for i in range(Ncholbeads):
            z_avg += y[i,2,t]

        z_avg = z_avg / N

        for i in range(Nlipidbeads):
            x[i,2,t] -= z_avg

        for i in range(Ncholbeads):
            y[i,2,t] -= z_avg

    return x,y
    
    
class Particle_d: #very specific to the displacement code
    def __init__(self,com):
        NDIM = len(com)
        Nconf = len(com[0])
        #nlin = Nconf/nlog
        self.pos = com
        #self.dr = np.zeros([NDIM,Nconf])
        self.s = np.zeros([Nconf])

    def calcS(self,t,frontblock,L):
        NDIM = len(self.pos)
        dr = [0,0,0]
        for k in range(NDIM):
            dr[k] = self.pos[k][t] - self.pos[k][frontblock]
            dr[k] = periodic(dr[k],L[k][t])
        self.s[t] = sqrt(dr[0]**2+dr[1]**2)
    
    def getS(self,t):
        return self.s[t]

def chunks(l, n): #specific to mobility g(r) code
    """Yield successive n-sized chunks from l."""
    """Small modification such that first and last chunk are of same size. Middle chunk is messed up"""
    size = len(l)//n
    lst = [size for i in range(n)]
    leftover = len(l) - size * n
    mid = n//2
    
    lst[mid] += leftover
    """
    if leftover < size//2:     
        lst[mid] += leftover
    else:
        lst.insert(mid+1,leftover)
    """
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
    
 
    """
for t in range(1,Nconf):
    if t%nlog == 0:
        frontblock = t
    for k in range(NDIM):
        dr = self.pos[k][t] - self.pos[k][frontblock]
        dr = periodic(dr,L[k][t])
        self.dr[k][t] = dr
    """     
            
def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError
