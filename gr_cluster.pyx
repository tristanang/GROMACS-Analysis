import lipid_cy as lib
import numpy as np
import sys

from math import sqrt
from math import pi

from jenkspy import jenks_breaks

def gr_cluster(trajFileName,Nconf,nlog,topology):
    #preprocessor
    NDIM = 3
    nlog = 46
    DR = 0.02

    cdef int i,j,k,z

    #cluster setting
    Nclusters = [2,3]
    times = [10,19,28,37,40,45]

    #Initializing the parameters
    Nblock = Nconf//nlog
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

    #calc displacements
    for start in range(Nblock):
        frontblock = start*nlog
        for t in times:
            endblock = frontblock + t
            for i in range(Nlipids):
                p_lipids[i].calcS(endblock,frontblock,L)
            for i in range(Nchol):
                p_chol[i].calcS(endblock,frontblock,L)
                
    #split upper, lower
    lipids = []
    lipids.append( list(filter(lambda x : x.pos[2][0] > 0 , p_lipids)) )
    lipids.append( list(filter(lambda x : x.pos[2][0] <= 0 , p_lipids)) )
    del p_lipids

    chol = []
    chol.append( list(filter(lambda x : x.pos[2][0] > 0 , p_chol)) )
    chol.append( list(filter(lambda x : x.pos[2][0] <= 0 , p_chol)) )
    del p_chol

    jenks_intervals = {}

    for start in range(Nblock):
        frontblock = start*nlog
        for t in times:
            endblock = frontblock + t
            jenks_intervals[endblock] = {}
           
            for clusterSize in Nclusters:
                jenks_intervals[endblock][clusterSize] = [0 for i in range(len(lipids))]
                
                for i in range(len(lipids)):    
                    jenks_intervals[endblock][clusterSize][i] = jenks_breaks(list(map(lambda x : x.getS(endblock),lipids[i])),clusterSize)

    for start in range(Nblock):
        frontblock = start*nlog
        for t in times:
            endblock = frontblock + t   

    #init gr dicts
    lipid_gr = {}
    cross_gr = {}
    lipid_pairs = {}
    cross_pairs = {}
    lipid_norm = {}
    cross_norm = {}

    for t in times:
        lipid_gr[t] = {}
        cross_gr[t] = {}
        lipid_pairs[t] = {}
        cross_pairs[t] = {}
        lipid_norm[t] = {}
        cross_norm[t] = {}
        
        for clusterSize in Nclusters:
            lipid_gr[t][clusterSize] = np.zeros((int(L[0][0]/DR)),np.int)
            cross_gr[t][clusterSize] = np.zeros((int(L[0][0]/DR)),np.int)
            lipid_pairs[t][clusterSize] = 0
            cross_pairs[t][clusterSize] = 0
    #

    def clusterForm(group,interval,endblock):
        clusters = []
        i_start = 0
        for i in range(1,len(interval)-1):
            i_end = lib.index(list(map(lambda x : x.getS(endblock),group)),interval[i])
            cluster = group[i_start:i_end+1]
            i_start = i_end + 1
            clusters.append(cluster)
        clusters.append(group[i_start:])
        
        return clusters
    #

    r = [0,0]
        
    for start in range(Nblock):
        frontblock = start*nlog
        for t in times:
            endblock = frontblock + t
        
            for z in range(len(lipids)): #upper lower
                lipids[z] = sorted(lipids[z], key = lambda x: x.getS(endblock))
                
                for clusterSize in Nclusters:
                    interval = jenks_intervals[endblock][clusterSize][z]
                    clusters = clusterForm(lipids[z],interval,endblock)
                    
                    for cluster in clusters:

                        
                        for i in range(len(cluster)):
                            for j in range(i+1,len(cluster)):
                                for k in range(NDIM-1):
                                    r[k] = cluster[i].pos[k][frontblock] - cluster[j].pos[k][frontblock]
                                    r[k] = lib.periodic(r[k],L[k][frontblock])
                                    
                                r2 = sqrt(r[0]*r[0]+r[1]*r[1])
                                lipid_gr[t][clusterSize][int(r2/DR)] += 1
                                lipid_pairs[t][clusterSize] += 1 
                                
                        for i in range(len(cluster)):
                            for j in range(len(chol[z])):
                                for k in range(NDIM-1):
                                    r[k] = cluster[i].pos[k][frontblock] - chol[z][j].pos[k][frontblock]
                                    r[k] = lib.periodic(r[k],L[k][frontblock])
                                    
                                r2 = sqrt(r[0]*r[0]+r[1]*r[1])
                                cross_gr[t][clusterSize][int(r2/DR)] += 1        
                                cross_pairs[t][clusterSize] += 1
    #printing

    L_ave = [np.mean(L[0]),np.mean(L[1])]

    for t in times:
        for clusterSize in Nclusters:
            
            lipid_norm[t][clusterSize] = lipid_pairs[t][clusterSize] * pi * DR/2.
            lipid_norm[t][clusterSize] = L_ave[0]*L_ave[1] / (4*lipid_norm[t][clusterSize])
            cross_norm[t][clusterSize] = cross_pairs[t][clusterSize] * pi * DR/2.
            cross_norm[t][clusterSize] = L_ave[0]*L_ave[1] / (4*cross_norm[t][clusterSize])

    size = (int(L[0][0]/DR))//2

    for t in times:
        for clusterSize in Nclusters:
            for speed in range(clusterSize):
            
                lipidFile = open("gr_m_lipids_c"+str(clusterSize)+"_t"+str(t)+"_s"+str(speed)+".dat",'w')
                for step in range(size):
                    lipidFile.write(str((step+0.5)*DR)+" "+str(lipid_norm[t][clusterSize]*lipid_gr[t][clusterSize][step]/((step+0.5)*DR))+"\n")
                
                lipidFile.close()
                
                cholFile = open("gr_m_cross_c"+str(clusterSize)+"_t"+str(t)+"_s"+str(speed)+".dat",'w')
                for step in range(size):
                    cholFile.write(str((step+0.5)*DR)+" "+str(cross_norm[t][clusterSize]*cross_gr[t][clusterSize][step]/((step+0.5)*DR))+"\n")
       
                cholFile.close()








