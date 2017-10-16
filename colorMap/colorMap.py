import sys
import re
import lipid as lib
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

#"preprocessor"
NDIM = 3

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
Nconf = 100
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

#Translating Z-Axis
x,y = lib.translateZ(x,y)

#COM MASSING
com_lipids = lib.comassing(x,"DPPC",L)
com_chol = lib.comassing(y,"CHOL",L)
del x,y

#Particle work
p_lipids = np.empty([Nlipids],dtype=object)
p_chol = np.empty([Nchol],dtype=object)

for i in range(Nlipids):
    p_lipids[i] = lib.Particle_d(com_lipids[i],L)

for i in range(Nchol):
    p_chol[i] = lib.Particle_d(com_chol[i],L)

#writing to file
for t in range(Nconf):
    if t%nlog == 0: #starting block check
        frontblock = t
        #technically rest of this if part is useless
        output_lipid = open("dr_lipid_top-"+str(frontblock)+"-xy.dat",'w')
        
        for i in range(Nlipids):
            if p_lipids[i].pos[2][frontblock] > 0:
                output_lipid.write(str(p_lipids[i].pos[0][frontblock])+" "+str(p_lipids[i].pos[1][frontblock])+"\n")
                
        output_lipid.close()

        output_chol = open("dr_chol_top-"+str(frontblock)+"-xy.dat",'w')
        
        for i in range(Nchol):
             if p_chol[i].pos[2][frontblock] > 0:
                output_chol.write(str(p_chol[i].pos[0][frontblock])+" "+str(p_chol[i].pos[1][frontblock])+"\n")

        output_chol.close()
        
    elif (t%nlog)%10 == 0 or t%nlog == (nlog - 1):
        for i in range(Nlipids):
            p_lipids[i].calcS(t,frontblock,L)
        
        for i in range(Nchol):
            p_chol[i].calcS(t,frontblock,L)

        outfile = open("dr_lipid_top-"+str(frontblock)+"-"+str(t)+".dat",'w')
        
        for i in range(Nlipids):
            if p_lipids[i].pos[2][frontblock] > 0:
                outfile.write(str(p_lipids[i].pos[0][frontblock])+" "+str(p_lipids[i].pos[1][frontblock])+" "+str(p_lipids[i].s[t])+"\n")
            
        outfile.close()

        outfile = open("dr_chol_top-"+str(frontblock)+"-"+str(t)+".dat",'w')
        
        for i in range(Nchol):
            if p_chol[i].pos[2][frontblock] > 0:
                outfile.write(str(p_chol[i].pos[0][frontblock])+" "+str(p_chol[i].pos[1][frontblock])+" "+str(p_chol[i].s[t])+"\n")
                
        outfile.close()
        """
        outfile = open("dr_lipid_top-"+str(t)+"-xy.dat",'w')
        
        for i in range(Nlipids):
            if p_lipids[i].pos[2][frontblock] > 0:
                outfile.write(str(p_lipids[i].pos[0][t])+" "+str(p_lipids[i].pos[1][t])+"\n")
        
        outfile.close()

        outfile = open("dr_chol_top-"+str(t)+"-xy.dat",'w')
        
        for i in range(Nchol):
            if p_chol[i].pos[2][frontblock] > 0:
                outfile.write(str(p_chol[i].pos[0][t])+" "+str(p_chol[i].pos[1][t])+"\n")

        outfile.close()
        """



def test(k,t):
    lst0 = []
    lst1 = []
    lst2 = []
    for i in range(Nlipids):
        lst0.append(p_lipids[i].pos[0][t])
        lst1.append(p_lipids[i].pos[1][t])           
        lst2.append(sqrt(p_lipids[i].dr[0][t]**2 + p_lipids[i].dr[1][t]**2))
    #plt.plot(lst2,'.')
    plt.scatter(lst0,lst1,c=lst2)
    plt.show()


