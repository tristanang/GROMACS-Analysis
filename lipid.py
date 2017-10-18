import re
import sys
import numpy as np
from math import sqrt

def cholConc(filename):
    file = open(filename,"r")
    for line in file:
        m = re.match("CHOL", line)
        if m:
            line = line.split()
            file.close()
            return int(line[1]) * 2

def processTraj(trajFileName,Nchol,NDIM,Nconf):
    trajFile = open(trajFileName,'r')
    N = int(trajFile.readline().split()[0])
    Nperlipid = 12
    Nperchol = 8
    Nlipids = (N - Nperchol * Nchol) // Nperlipid
    Ncholbeads = Nchol * Nperchol
    Nlipidbeads = Nlipids * Nperlipid
    
    #reading traj file
    trajFile.close()
    trajFile = open(trajFileName,'r')

    L = np.empty([NDIM,Nconf])
    x = np.empty([Nlipidbeads,NDIM,Nconf])
    y = np.empty([Ncholbeads,NDIM,Nconf])

    for t in range(Nconf):
        trajFile.readline()
        
        #Box Sizes
        for k in range(NDIM):
            L[k][t] = trajFile.readline().strip()
        
        #Bead Coordinates
        for i in range(Nlipidbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
                
            for k in range(NDIM):
                x[i][k][t] = config[k]
            
        for i in range(Ncholbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
            
            for k in range(NDIM):
                y[i][k][t] = config[k]
            
        for i in range(Nlipidbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
             
            for k in range(NDIM):
                x[i+Nlipidbeads//2][k][t] = config[k]
        
        for i in range(Ncholbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
            
            for k in range(NDIM):
                y[i+Ncholbeads//2][k][t] = config[k]

    trajFile.close()

    return N,L,x,y
    
def comassing(ar,bead,L=None):
    if bead == "DPPC":
        size = 12
        #center = 2
        #3rd bead but zero indexed
    elif bead == "CHOL":
        size = 8
        #center = 4
    else:
        print("bead type error")

    N = len(ar)//size
    NDIM = len(ar[0])
    Nconf = len(ar[0][0])

    comAr = np.zeros([N,NDIM,Nconf])

    for t in range(Nconf):
        for k in range(NDIM):
            for i in range(N):
                ii = i*size
                
                #comAr[i][k][t] = ar[ii+center][k][t]
                
                for j in range(size): #where's the 0 bead
                    ar[ii+j][k][t] = periodic_p(ar[ii+j][k][t],ar[ii][k][t],L[k][t])
                    #fix your boundary condition
                    #look at structure factor
                    comAr[i][k][t] += ar[ii+j][k][t]/size
  
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

def translateZ(x,y):
    Nconf = len(x[0][0])
    Nlipidbeads = len(x)
    Ncholbeads = len(y)
    N = Nlipidbeads + Ncholbeads

    for t in range(Nconf):
        z_avg = 0.
            
        for i in range(Nlipidbeads):
            z_avg += x[i][2][t]

        for i in range(Ncholbeads):
            z_avg += y[i][2][t]

        z_avg = z_avg / N

        for i in range(Nlipidbeads):
            x[i][2][t] -= z_avg

        for i in range(Ncholbeads):
            y[i][2][t] -= z_avg

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
            
