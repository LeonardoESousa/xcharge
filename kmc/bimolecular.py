import numpy as np
import random
from kmc.particles import *


# BIMOLEC FUNCS NOTE: ALL THEM MUST HAVE THE SAME VARIABLES (system,tipos,Ss,indices,locs)

#recombination electron-hole pair
def ele_hol_recomb(Ss,system,superp):
    if random.uniform(0,1) <= 0.75 and abs(Ss[0].identity) != abs(Ss[1].identity):
        system.set_particles([Triplet(Ss[0].position)])
    else:
        system.set_particles([Singlet(Ss[0].position)])

    for i in superp:
        Ss[i].kill('recomb',system,system.lumo)

           
#singlet-singlet annihilation (ssa)   
def anni_sing(Ss,system,superp):
    Ss[random.choices(superp)[0]].kill('ssa',system,system.s1)
    

bimolec_funcs_array = {('singlet','singlet'):anni_sing, ('electron','hole'):ele_hol_recomb}   
############################################           