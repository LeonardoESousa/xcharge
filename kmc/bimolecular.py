import numpy as np
import random
from kmc.particles import *


# BIMOLEC FUNCS NOTE: ALL THEM MUST HAVE THE SAME VARIABLES (system,tipos,Ss,indices,locs)

#recombination electron-hole pair
def ele_hol_recomb(Ss,system,superp):
    first, second = (Ss[int(index)] for index in superp)
    if random.uniform(0,1) <= 0.75 and abs(first.identity) != abs(second.identity):
        system.set_particles([Triplet(first.position)])
    else:
        system.set_particles([Singlet(first.position)])

    for i in superp:
        Ss[i].kill('recomb',system,system.lumo,'converted')

           
#singlet-singlet annihilation (ssa)   
def anni_sing(Ss,system,superp):
    Ss[random.choices(superp)[0]].kill('ssa',system,system.s1,'dead')


def formation(Ss,system,superp):
    first = Ss[int(superp[0])]
    system.set_particles([Frenkelpair(first.position)])

    for i in superp:
        Ss[i].kill('recomb',system,system.lumo,'converted')

bimolec_funcs_array = {('singlet','singlet'):anni_sing, ('electron','hole'):ele_hol_recomb, ('interstitial','vacancy'):formation, ('vacancy','interstitial'):formation}   
############################################           
