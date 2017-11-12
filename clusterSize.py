import lipid_cy as lib
import numpy as np
from jenkspy import jenks_breaks

from bisect import bisect_left

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError
    
def jenksSize(p_lipids,clustersize,t):
    lipids = list(map(lambda x : x.getS(t),p_lipids))
    lipids.sort()
    interval = jenks_breaks(lipids,clustersize)
    lst = []
    i_start = 0
    for i in range(1,len(interval)):
        i_end = index(lipids,interval[i])
        lst.append((i_end-i_start+1)/len(lipids))
        i_start = i_end + 1
    
    return lst

#preprocessor
NDIM = 3
nlog = 46

#user input
trajFileName = "testData/temp20"
Nconf = 92
topology = "testData/20chol.top"

#cluster setting
Nclusters = [2,3,4,5]

#Initializing the parameters
Nblock = Nconf//nlog
Nchol = lib.cholConc(topology)
N,L,x,y = lib.processTraj(trajFileName,Nchol,NDIM,Nconf)

print("hi")

Nperlipid = 12
Nperchol = 8
Nlipids = (N - Nperchol * Nchol) // Nperlipid
Ncholbeads = Nchol * Nperchol
Nlipidbeads = Nlipids * Nperlipid

#Translating z-axis
x,y = lib.translateZ(x,y)

print("hi")

#COM MASSING
com_lipids = lib.comassing(x,"DPPC",L)
com_chol = lib.comassing(y,"CHOL",L)
del x,y

print("hi")

#Particle Class
p_lipids = [0 for i in range(Nlipids)]
p_chol = [0 for i in range(Nchol)]

for i in range(Nlipids):
    p_lipids[i] = lib.Particle_d(com_lipids[i])
for i in range(Nchol):
    p_chol[i] = lib.Particle_d(com_chol[i])

del com_lipids, com_chol

#calc displacements
for t in range(Nconf):
    if t%46 == 0:
        frontblock = t
    else:
        for i in range(Nlipids):
            p_lipids[i].calcS(t,frontblock,L)

#split upper, lower
lipids_upper = list(filter(lambda x : x.pos[2][0] > 0 , p_lipids)) 
lipids_lower = list(filter(lambda x : x.pos[2][0] <= 0 , p_lipids)) 
del p_lipids
Nlipids_upper = len(lipids_upper) 
Nlipids_lower = len(lipids_lower)

#chol_upper = 
#p_chol = list(filter(lambda x : x.pos[2][0] > 0 , p_chol))
#Nlipids = len(p_lipids)
#Nchol = len(p_chol)
#N = Nlipids + Nchol

#intervals
output = {}

for clustersize in Nclusters:
    output[clustersize] = np.zeros([clustersize,nlog])

for clustersize in Nclusters:
    for t in range(Nconf):
        if t%nlog == 0:
            frontblock = t 
        else:
            lst1 = jenksSize(lipids_upper,clustersize,t)
            lst2 = jenksSize(lipids_lower,clustersize,t)
            for i in range(clustersize):
                output[clustersize][i][t%nlog] += lst1[i]/(2.*Nblock)
                output[clustersize][i][t%nlog] += lst2[i]/(2.*Nblock)

#printing
                
for clustersize in Nclusters:
    outFile = open("clusterSize="+str(clustersize)+".dat","w")
    for t in range(1,nlog): #should i do 1,nlog to remove 0?
        write = ""
        for i in range(clustersize):
            write += str(output[clustersize][i][t]) + " "
        outFile.write(str(write)+"\n")
   
        




