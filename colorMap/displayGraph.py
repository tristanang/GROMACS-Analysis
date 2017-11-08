import sys
import re
import lipid as lib
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

def display(dr_l,lipid_xy,chol_xy,dr_c=""):
    chol_x = []
    chol_y = []
    lipid_x = []
    lipid_y = []
    dr_lipid = []
    dr_chol = []

    file = open(dr_l,'r')
    for line in file:
        c1,c2,c3 = map(float, line.split())
        dr_lipid.append(c3)
    file.close()

    file = open(lipid_xy,'r')
    for line in file:
        c1,c2 = map(float,line.split())
        lipid_x.append(c1)
        lipid_y.append(c2)
    file.close()

    file = open(chol_xy,'r')
    for line in file:
        c1,c2 = map(float,line.split())
        chol_x.append(c1)
        chol_y.append(c2)
    file.close()

    plt.scatter(lipid_x,lipid_y,c=dr_lipid,cmap='viridis')
    plt.colorbar()
    plt.scatter(chol_x,chol_y,color='black')
    plt.show()


display("dr_lipid_top-0-45.dat","dr_lipid_top-0-xy.dat","dr_chol_top-0-xy.dat")
