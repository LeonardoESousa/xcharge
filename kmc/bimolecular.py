import random
import kmc.particles


# BIMOLEC FUNCS NOTE: ALL THEM MUST HAVE THE SAME VARIABLES (system,tipos,Ss,indices,locs)


# recombination electron-hole pair
def ele_hol_recomb(Ss, system, superp):
    if random.uniform(0, 1) <= 0.75 and abs(Ss[0].identity) != abs(Ss[1].identity):
        system.set_particles([kmc.particles.Triplet(Ss[0].position)])
    else:
        system.set_particles([kmc.particles.Singlet(Ss[0].position)])

    for i in superp:
        Ss[i].kill("recomb", system, system.lumo, "converted")


# singlet-singlet annihilation (ssa)
def anni_sing(Ss, system, superp):
    s = Ss[random.choices(superp)[0]]
    s.kill("ssa", system, 0, "dead")


bimolec_funcs_array = {
    ("singlet", "singlet"): anni_sing,
    ("electron", "hole"): ele_hol_recomb,
}
############################################
