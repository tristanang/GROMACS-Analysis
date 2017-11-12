import sys
import re
import lipid_cy as lib
import numpy as np
from math import sqrt
from math import pi
import matplotlib.pyplot as plt
from bisect import bisect_left

#"preprocessor"
		
NDIM = 3
DR = 0.02
dt = 0.02
confGroups = 5
Nchunks = 5 #with this number
"""
#Arguments of script
if len(sys.argv) != 5:
    print("You need five arguments")
    print("Trajectory File")
    print("Number of Configurations")
    print("Logarithmic Block Size")
    print("Cholesterol Topology File")
    exit()

trajFileName = sys.argv[1]
Nconf = int(sys.argv[2])
nlog = int(sys.argv[3])
topology = sys.argv[4]

"""
trajFileName = "temp20"
Nconf = 92
nlog = 46
topology = "20chol.top"


#Initializing the parameters
Nchol = lib.cholConc(topology)
N,L,x,y = lib.processTraj(trajFileName,Nchol,NDIM,Nconf)

Nperlipid = 12
Nperchol = 8
Nlipids = (N - Nperchol * Nchol) // Nperlipid
Ncholbeads = Nchol * Nperchol
Nlipidbeads = Nlipids * Nperlipid

#Translating z-axis
x,y = lib.translateZ(x,y)

#COM MASSING
com_lipids = lib.comassing(x,"DPPC",L)
com_chol = lib.comassing(y,"CHOL",L)
del x,y

#Particle Class
p_lipids = [0 for i in range(Nlipids)]
p_chol = [0 for i in range(Nchol)]

for i in range(Nlipids):
    p_lipids[i] = lib.Particle_d(com_lipids[i])
for i in range(Nchol):
    p_chol[i] = lib.Particle_d(com_chol[i])

del com_lipids, com_chol 
 
for t in range(Nconf):
    if t%nlog == 0: #starting block check
        frontblock = t
    elif (t%nlog)%10 == 0 or t%nlog == (nlog - 1):
            for i in range(Nlipids):
                p_lipids[i].calcS(t,frontblock,L)
            #for i in range(Nchol):
                #p_chol[i].calcS(t,frontblock,L)

lipid_gr = np.zeros([confGroups,Nchunks,int(L[0][0]/DR)],int)
cross_gr =  np.zeros([confGroups,Nchunks,int(L[0][0]/DR)],int)
r = [0,0]
lipid_pair = np.zeros([confGroups,Nchunks])
cross_pair = np.zeros([confGroups,Nchunks])
                
for t in range(Nconf):
    if t%nlog == 0:
        frontblock = t
        confCount = 0
    elif (t%nlog)%10 == 0 or t%nlog == (nlog - 1):
        p_lipids = sorted(p_lipids, key = lambda x: x.getS(t))
        #p_chol = sorted()
        
        chunkGen = lib.chunks(p_lipids,Nchunks)
        chunkCount = 0
        for chunk in chunkGen:
            
            for i in range(len(chunk)):
                for j in range(i+1,len(chunk)):
                    if chunk[i].pos[2][frontblock]*chunk[j].pos[2][frontblock] > 0: 
                        for k in range(NDIM-1):
                            r[k] = chunk[i].pos[k][frontblock] - chunk[j].pos[k][frontblock]
                            r[k] = lib.periodic(r[k],L[k][frontblock])
                            
                        r2 = sqrt(r[0]*r[0] + r[1]*r[1])
                        lipid_gr[confCount][chunkCount][int(r2/DR)] += 1
                        lipid_pair[confCount][chunkCount] += 1
                        
            for i in range(len(chunk)):
                for j in range(Nchol):
                    if chunk[i].pos[2][frontblock]*p_chol[j].pos[2][frontblock] > 0:
                        for k in range(NDIM-1):
                            r[k] = chunk[i].pos[k][frontblock] - p_chol[j].pos[k][frontblock]
                            r[k] = lib.periodic(r[k],L[k][frontblock])
                            
                        r2 = sqrt(r[0]*r[0] + r[1]*r[1])
                        #if r2/DR < 1:
                         #   print(i,j,frontblock)
                        cross_gr[confCount][chunkCount][int(r2/DR)] += 1
                        cross_pair[confCount][chunkCount] += 1
                
            chunkCount += 1
            
        confCount += 1
        confCount %= confGroups

#printing
lipid_norm = np.zeros([confGroups,Nchunks])
cross_norm = np.zeros([confGroups,Nchunks])
L_ave = [np.mean(L[0]),np.mean(L[1])]

for i in range(confGroups):
    for j in range(Nchunks):
        lipid_norm[i][j] = lipid_pair[i][j] * pi * DR / 2.
        lipid_norm[i][j] = L_ave[0]*L_ave[1]/(4*lipid_norm[i][j])
        cross_norm[i][j] = cross_pair[i][j] * pi * DR /2.
        cross_norm[i][j] = L_ave[0]*L_ave[1]/(4*cross_norm[i][j])

for i in range(confGroups):
    for j in range(Nchunks):
        lipidFile = open("gr_m_lipids_t="+str(i)+"_speed="+str(j)+".dat",'w')
        for step in range(len(lipid_gr[i][j])//2):
            lipidFile.write(str((step+0.5)*DR)+" "+str(lipid_norm[i][j]*lipid_gr[i][j][step]/((step+0.5)*DR))+"\n")
        lipidFile.close()
        
        cholFile = open("gr_m_cross_t="+str(i)+"_speed="+str(j)+".dat",'w')
        for step in range(len(cross_gr[i][j])//2):
            cholFile.write(str((step+0.5)*DR)+" "+str(cross_norm[i][j]*cross_gr[i][j][step]/((step+0.5)*DR))+"\n")       
        cholFile.close()
        

                
                
                
                
                
                
