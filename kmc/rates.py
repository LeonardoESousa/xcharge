
import numpy as np
import random
from kmc.particles import *
import kmc.utils

epsilon_vaccum = 8.854187e-12    #Permitivity in C/Vm
e              = -1.60217662e-19 #Electron charge    
kb             = 8.617e-5        #Boltzmann constant
hbar           = 6.582e-16       #Reduced Planck's constant



###RATES#################################################################################    

##FUNCTION FOR SETTING RADII#############################################################

def raios(num,Rf, mat, lifetime, mats):
  # Initialize the Raios array with the value of Rf[(mat,mat)]
  Raios = np.empty(num)
  Raios.fill(Rf[(mat,mat)])

  # Use NumPy's where function to set the values of Raios for the other materials
  for m in lifetime.keys():
    if m != mat:
      Raios = np.where(mats == m, Rf[(mat,m)], Raios)

  return Raios


def raios_dist(num,Rf,mat,lifetime,mats):
    Raios = np.array(random.choices(Rf[(mat,mat)][:,0],Rf[(mat,mat)][:,1],k=num))
    materiais = [i for i in lifetime.keys() if i != mat]
    for m in materiais:
        R2 = np.array(random.choices(Rf[(mat,m)][:,0],Rf[(mat,m)][:,1],k=num))
        Raios[mats == m] = R2[mats == m]
    return Raios


#function to convert dictionary with (i,j) keys to ixj array
def dict_to_array(d):
    keys = d.keys()
    num_keys = len(set(key[0] for key in keys))
    radius = np.empty((num_keys,num_keys))
    for key in keys:
        radius[key[0],key[1]] = d[key]
    return radius

#########################################################################################

##STANDARD FORSTER TRANSFER RATE#########################################################
class Forster:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = dict_to_array(kwargs['Rf'])
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53


    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = kwargs['mats']
        local  = ex.position  
        mat    = kwargs['matlocal']
        num = len(mats)
        taxa = kmc.utils.forster(self.Rf[mat,:],mats,num,self.alpha*self.mu[mat], r,1/self.lifetime[mat])
        return taxa


    def action(self,particle,system,local):
        particle.move(local,system)
#########################################################################################

##TRIPLET TO SINGLET FORSTER TRANSFER####################################################
class ForsterT:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = dict_to_array(kwargs['Rf'])
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53

    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = kwargs['mats']
        local  = ex.position  
        mat    = kwargs['matlocal']
        num = len(mats)
        taxa = kmc.utils.forster(self.Rf[mat,:],mats,num,self.alpha*self.mu[mat], r,1/self.lifetime[mat])
        return taxa

    def action(self,particle,system,local):
        particle.move(local,system)
        particle.kill('tts',system,system.t1,'converted')
        system.set_particles([Singlet(local)])
#########################################################################################

##FORSTER ANNIHILATION RADIUS#########################################################
class Forster_Annirad:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rf = dict_to_array(kwargs['Rf'])
        self.lifetime = kwargs['life']
        self.mu = kwargs['mu']
        self.alpha = 1.15*0.53
        self.anni_rad = kwargs['anni_rad']

    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = kwargs['mats']
        cut    = kwargs['cut']
        local  = ex.position  
        mat    = kwargs['matlocal']
        num = len(mats)
        relevant_particles = [p for p in system.particles if p.identity != ex.identity and p.position in cut]
        ss = [(np.where(cut == p.position)[0][0],self.anni_rad[(mat,system.mats[p.position])][p.species]) for p in relevant_particles]
        replace_pos   = np.array([ele[0] for ele in ss],dtype=np.int32)
        replace_raios = np.array([ele[1] for ele in ss],dtype=np.double)
        mum = len(replace_pos)
        taxa = kmc.utils.forster_anni(self.Rf[mat,:],mats,num,self.alpha*self.mu[mat], r,1/self.lifetime[mat], replace_pos, replace_raios, mum)
        return taxa


    def action(self,particle,system,local):
        particle.move(local,system)       
#########################################################################################

##STANDARD DEXTER TRANSFER RATE##########################################################
class Dexter:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.Rd = dict_to_array(kwargs['Rd'])
        self.lifetime = kwargs['life']
        self.L = kwargs['L']
       
    def rate(self,**kwargs):
        r      = kwargs['r']
        system = kwargs['system']
        ex     = kwargs['particle']
        mats   = kwargs['mats']
        local  = ex.position  
        mat    = kwargs['matlocal']
        num = len(mats)
        taxa = kmc.utils.dexter(self.Rd[mat,:],1/self.L[mat],1/self.lifetime[mat], mats,num, r)
        return taxa

    def action(self,particle,system,local):
        particle.move(local,system)
#########################################################################################       

##EXCITON DISSOCIATION BY ELECTRON HOP RATE##############################################
class Dissociation_electron:
    def __init__(self,**kwargs):
        self.kind = 'dissociation_e'
        self.AtH = kwargs['AtH']
        self.inv      = kwargs['invrad']
        self.T        = kwargs['T']

    def rate(self,**kwargs):
        system   = kwargs['system']
        r        = kwargs['r']
        particle = kwargs['particle']
        local    = particle.position
        cut      = kwargs['cut']
        mats   = kwargs['mats']
        mat    = kwargs['matlocal']
        num      = len(mats)

        lumos = system.lumo
        homos = system.homo
        if particle.species   == 'singlet':
            s1s   = system.s1 
        elif particle.species == 'triplet':
            s1s   = system.t1

        AtH        = raios(num,self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]
        
        DEe = lumos[cut] - (homos[local] + s1s[local])
        taxae = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r -(DEe+abs(DEe))/(2*kb*self.T))
        taxae[r == 0] = 0
        return taxae
            
     
    def action(self,particle,system,local):
        e = Electron(local)
        h = Hole(particle.position)
        h.identity = -1*e.identity   
        system.set_particles([e,h])
        particle.kill(self.kind,system,system.s1,'converted')
#########################################################################################

##EXCITON DISSOCIATION BY HOLE HOP RATE##################################################
class Dissociation_hole:
    def __init__(self,**kwargs):
        self.kind = 'dissociation_h'
        self.AtH = kwargs['AtH']
        self.inv      = kwargs['invrad']
        self.T        = kwargs['T']

    def rate(self,**kwargs):
        system   = kwargs['system']
        r        = kwargs['r']
        particle = kwargs['particle']
        local    = particle.position 
        cut      = kwargs['cut']
        mats   = kwargs['mats']
        mat    = kwargs['matlocal']
        num      = len(mats)
        
        lumos = system.lumo
        homos = system.homo
        if particle.species   == 'singlet':
            s1s   = system.s1 
        elif particle.species == 'triplet':
            s1s   = system.t1

        AtH        = raios(num,self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]
        
        DEh = (lumos[local] - s1s[local]) - homos[cut]
        taxah = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r -(DEh+abs(DEh))/(2*kb*self.T))
        taxah[r == 0] = 0
        return taxah
                 
    def action(self,particle,system,local):
        e = Electron(particle.position)
        h = Hole(local) 
        h.identity = -1*e.identity   
        system.set_particles([e,h])
        particle.kill(self.kind,system,system.s1,'converted')
#########################################################################################


def corrected_energies(system,s,r,dx,dy,dz,cut):
    electrostatic = system.electrostatic()
    potential = electrostatic[cut] - (
        system.field[0] * dx + system.field[1] * dy + system.field[2] * dz
    )
    r_safe = np.where(r != 0, r, np.inf)
    potential -= s.charge * abs(e) / (4 * np.pi * system.epsilon * r_safe)
    electron_sites = set(system.species_positions("electron")) - {s.position}
    hole_sites = set(system.species_positions("hole")) - {s.position}
    electron_mask = np.fromiter(
        (int(site) in electron_sites for site in cut), dtype=bool, count=len(cut)
    )
    hole_mask = np.fromiter(
        (int(site) in hole_sites for site in cut), dtype=bool, count=len(cut)
    )
    origin_potential = electrostatic[s.position]

    if s.species == 'electron':
        engs = np.array(system.lumo[cut], copy=True)
        engs[hole_mask] = system.homo[cut[hole_mask]]
        engs[electron_mask] = np.inf
        engs -= potential
        origin_energy = system.lumo[s.position] - origin_potential
        difference = engs - origin_energy
        DE = difference + abs(difference)
    elif s.species == 'hole':
        engs = np.array(system.homo[cut], copy=True)
        engs[electron_mask] = system.lumo[cut[electron_mask]]
        engs[hole_mask] = -np.inf
        engs -= potential
        origin_energy = system.homo[s.position] - origin_potential
        difference = origin_energy - engs
        DE = difference + abs(difference)
    return DE



##MILLER-ABRHAMS RATE####################################################################
class MillerAbrahams:
    def __init__(self,**kwargs):
        self.kind = 'miller-abrahams'
        self.AtH  = kwargs['AtH']
        self.inv  = kwargs['invrad']
        self.T    = kwargs['T']

    def rate(self,**kwargs):
        system    = kwargs['system']
        r         = kwargs['r']
        dx        = kwargs['dx']
        dy        = kwargs['dy']
        dz        = kwargs['dz']
        particle  = kwargs['particle']
        cut       = kwargs['cut']
        mats   = kwargs['mats']
        mat    = kwargs['matlocal']
        num      = len(mats)        
        
        AtH        = raios(len(r),self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]

        DE = corrected_energies(system,particle,r,dx,dy,dz,cut)
        taxa = (1e-12)*(AtH)*np.exp(-2*in_loc_rad*r -DE/(2*kb*self.T)) 
        
        taxa[r == 0] = 0
        return taxa
 
    def action(self,particle,system,local):
        particle.move(local,system)

######################################################################################### 



#MONOMOLECULAR RATES#####################################################################

##FLUORESCENCE RATE######################################################################
class Fluor:
    def __init__(self,**kwargs):
        self.kind = 'fluor'
        self.lifetime = kwargs['life']

    def rate(self,**kwargs):
        return 1/self.lifetime[kwargs['material']]
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.s1,'dead') 
#########################################################################################

##PHOSPHORESCENCE RATE###################################################################
class Phosph:
    def __init__(self,**kwargs):
        self.kind = 'phosph'
        self.lifetime = kwargs['life']

    def rate(self,**kwargs):
        return 1/self.lifetime[kwargs['material']]
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.t1,'dead')
#########################################################################################
 
##NONRADIATIVE DECAY RATE################################################################         
class Nonrad:
    def __init__(self,**kwargs):
        self.kind = 'nonrad'
        self.taxa = kwargs['rate']

    def rate(self,**kwargs):
        return self.taxa[kwargs['material']]
     
    def action(self,particle,system,local):
        particle.kill(self.kind,system,system.s1,'dead')        
#########################################################################################

##INTERSYSTEM CROSSING RATE##############################################################
class ISC:
    def __init__(self,**kwargs):
        self.kind = 'isc'
        self.taxa = kwargs['rate']
        self.map  = {'singlet':'triplet', 'triplet':'singlet'}

    def rate(self,**kwargs):
        material = kwargs['material']
        return self.taxa[material]
     
    def action(self,particle,system,local):
        if particle.species == 'singlet':
            system.set_particles([Triplet(particle.position)])
            particle.kill(self.kind,system,system.s1,'converted')
        elif particle.species == 'triplet':    
            system.set_particles([Singlet(particle.position)])
            particle.kill('r'+self.kind,system,system.s1,'converted')
#########################################################################################


def filter(num,ks, mat, mats, materials_list):
    rates = np.zeros(num, dtype=float)
    for material in materials_list:
        if material == 999:
            continue
        rates[mats == material] = ks.get((mat, material), 0.0)
    return rates


##DEFECTS MIGRATION RATE##############################################################
class Migration:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.neighbor_mode = 'first'
        self.k = dict(kwargs['k'])
        self.materials_list = list({key[1] for key in self.k})

    def rate(self,**kwargs):
        system   = kwargs['system']
        mats     = kwargs['mats']    
        mat      = kwargs['matlocal']
        cut      = kwargs['cut']   # full-lattice indices of those neighbors

        #migration can not occur in occupied sites!
        occupied = np.fromiter(
            (system.occupancy.get(int(site), 0) > 0 for site in cut),
            dtype=bool,
            count=len(cut),
        )
        taxa = filter(len(mats),self.k,mat,mats,self.materials_list)
        return np.where(occupied, 0.0, taxa)


    def action(self,particle,system,local):
        particle.move(local,system)
class DissociationFP:
    def __init__(self,**kwargs):
        self.kind = 'dissociation_fp'
        self.k = dict(kwargs['k'])
        self.k[999] = 0
    def rate(self,**kwargs):
        return self.k[kwargs['material']]
                 
    def action(self,particle,system,local):
        #convention:
        # ->local = particle.position where the FP is centered, conventioned to be the place where the Intersitial moves to
        # ->ghost  position is reserved for the vacancy
        I = Interstitial(particle.position)
        V = Vacancy(particle.ghost_site)
        system.mats[particle.ghost_site] = particle.origin_site 
        system.set_particles([I,V])
        #print(f'diss:vac->{particle.ghost_site} inter->{particle.position}')
        particle.kill(self.kind,system,system.s1,'converted')
#########################################################################################        

class Annihilation:
    def __init__(self,**kwargs):
        self.kind = 'annihilation'
        self.k = dict(kwargs['k'])
        self.k[999] = 0
    def rate(self,**kwargs):
        return self.k[kwargs['material']]
     
    def action(self,particle,system,local):
        #print('ani:',particle.species,particle.ghost_site,particle.origin_site,particle.status)
        system.mats[particle.ghost_site] = particle.origin_site #reverting the ghost site
        particle.kill(self.kind,system,system.s1,'dead') 
#########################################################################################
class Formation:
    def __init__(self, **kwargs):
        self.kind = 'formation'
        self.neighbor_mode = 'first'
        supplied_rate = kwargs['k']
        if np.isscalar(supplied_rate):
            self.k = None
            self.scalar_rate = float(supplied_rate)
            self.materials_list = []
        else:
            self.k = dict(supplied_rate)
            self.scalar_rate = None
            self.materials_list = list({key[1] for key in self.k})
    def rate(self, **kwargs):
        system   = kwargs['system']
        particle = kwargs['particle']
        mats     = kwargs['mats']    
        mat      = kwargs['matlocal']
        cut      = kwargs['cut']   # full-lattice indices of those neighbors

        origin_type = particle.species
        if origin_type == 'interstitial':
            neigh_type = 'vacancy'
        elif origin_type == 'vacancy':
            neigh_type = 'interstitial'        
        else:
            return np.zeros(len(cut), dtype=float)

        neighbor_positions = system.species_positions(neigh_type)
        occupied = np.fromiter(
            (int(site) in neighbor_positions for site in cut),
            dtype=bool,
            count=len(cut),
        )

        # Rate: k when destination has an insterstitial, else 0
        if self.k is None:
            taxa = np.full(len(mats), self.scalar_rate, dtype=float)
            taxa[mats == 999] = 0.0
        else:
            taxa = filter(len(mats),self.k,mat,mats,self.materials_list)
        taxa = np.where(occupied, taxa, 0.0)

        return taxa


    def action(self, particle, system, local):
        # FP generation requires the addition of a ghost(deactivated) site
        # where new particles cant hop to it 
        # FP.ghost_site is the position (in len(system.X))
        # material type 999 is reserved for this special site
        #convention:
        # ->local where the FP is centered, conventioned to be the place where the Intersitial moves to
        # ->ghost position is reserved for the vacancy
        
        origin_type = particle.species
        #two forms of formation: I first  or V first
        if origin_type == 'interstitial':
            neigh_type = 'vacancy'
        if origin_type == 'vacancy':
            neigh_type = 'interstitial'         
             
        if origin_type == 'interstitial':
          FP_position   = particle.position
          Ghost_position = local
        if origin_type == 'vacancy':
          FP_position    = local
          Ghost_position = particle.position
        
        FP = Frenkelpair0(FP_position)
        FP.ghost_site = Ghost_position
        FP.origin_site = system.mats[Ghost_position]
        system.mats[Ghost_position] = 999
        
        
        system.set_particles([FP])
        particle.kill(self.kind, system, system.s1, 'converted')
        #debug
        #xfp,yfp,zfp = system.X[particle.position],system.Y[particle.position],system.Z[particle.position]
        #xg,yg,zg = system.X[local],system.Y[local],system.Z[local]
        #dx,dy,dz =  xg-xfp, yg-yfp, zg-zfp
        #dist = np.sqrt(dx**2 + dy**2+dz**2)
        #print('form:',dist,particle.position,local)
        inters = [ p for p in system.particles if p.species==neigh_type]
        for p in inters:
            if p.position == local:
                p.kill(self.kind, system, system.s1, 'converted')
        '''
        for p in system.particles:
            if isinstance(p, Interstitial) and p.position == local:
                p.kill('interstitial', system, system.s1, 'converted')
                break
        '''
#########################################################################################
# I + I -> I2
class I2Formation:
    def __init__(self, **kwargs):
        self.kind = 'formation'
        self.neighbor_mode = 'first'
        self.k = dict(kwargs['k'])
        self.materials_list = list({key[1] for key in self.k})
    def rate(self, **kwargs):
        system   = kwargs['system']
        mats     = kwargs['mats']    
        mat      = kwargs['matlocal']
        cut      = kwargs['cut']   # full-lattice indices of those neighbors

        interstitial_sites = system.species_positions('interstitial')
        occupied = np.fromiter(
            (int(site) in interstitial_sites for site in cut),
            dtype=bool,
            count=len(cut),
        )

        # Rate: k when destination has an insterstitial, else 0
        taxa = filter(len(mats),self.k,mat,mats,self.materials_list)
        taxa = np.where(occupied, taxa, 0.0)

        return taxa


    def action(self, particle, system, local):
        I2_part = I2(local)
        I2_part.ghost_site = particle.position
        I2_part.origin_site = system.mats[particle.position]
        system.mats[particle.position] = 999
        system.set_particles([I2_part])
        particle.kill(self.kind, system, system.s1, 'converted')
        inters = [ p for p in system.particles if p.species=="interstitial"]
        for p in inters:
            if p.position == local:
                p.kill(self.kind, system, system.s1, 'converted')

class I2Dissociation:
    def __init__(self,**kwargs):
        self.kind = 'dissociation_I2'
        self.k = dict(kwargs['k'])
        self.k[999] = 0
    def rate(self,**kwargs):
        return self.k[kwargs['material']]
                 
    def action(self,particle,system,local):
        #convention:
        # ->local = particle.position where the I2 is centered, conventioned to be the place where the Intersitial moves to
        # ->ghost  position is reserved for the other I
        IA = Interstitial(particle.position)
        IB = Interstitial(particle.ghost_site)
        system.mats[particle.ghost_site] = particle.origin_site 
        system.set_particles([IA,IB])
        #print(f'diss:vac->{particle.ghost_site} inter->{particle.position}')
        particle.kill(self.kind,system,system.s1,'converted')    
####################################################################################
class FP_generation:
    def __init__(self,**kwargs):
        self.kind = 'generation'
        self.k = kwargs['k']
        self.pairs = kwargs['pairs']
        self.causamortis = kwargs.get('causamortis', "generated")

    def available_pairs(self, system):
        occupied = system.occupancy
        available = []
        for vacancy_material, interstitial_material in self.pairs:
            interstitial_sites = np.flatnonzero(
                system.mats == interstitial_material
            )
            for interstitial_site in interstitial_sites:
                interstitial_site = int(interstitial_site)
                if occupied.get(interstitial_site, 0):
                    continue
                vacancy_sites = system.first_neighbors_of_material(
                    interstitial_site, vacancy_material
                )
                for vacancy_site in vacancy_sites:
                    vacancy_site = int(vacancy_site)
                    if not occupied.get(vacancy_site, 0):
                        available.append((interstitial_site, vacancy_site))
        return available

    def rate(self, system=None, **kwargs):
        if self.k <= 0:
            return 0.0
        if system is not None and not self.available_pairs(system):
            return 0.0
        return float(self.k)
    def create(self,system,selected,ghost,**kwargs):
        FP = Frenkelpair2(selected)
        FP.ghost_site = ghost
        FP.origin_site = system.mats[ghost] #storing the old material type
        system.mats[ghost] = 999 #changing the mat
        system.set_particles([FP])
        causamortis = self.causamortis
        FP.make_text(system,system.s1,causamortis=causamortis)
        FP.stamp_time(system)
        
        #debug
        #xfp,yfp,zfp = system.X[selected],system.Y[selected],system.Z[selected]
        #xg,yg,zg = system.X[ghost],system.Y[ghost],system.Z[ghost]
        #dx,dy,dz =  xg-xfp, yg-yfp, zg-zfp
        #dist = np.sqrt(dx**2 + dy**2+dz**2)
        #print(self.causamortis,dist)
    def select(self,system,**kwargs):
        num = kwargs['num']
        candidates = self.available_pairs(system)
        random.shuffle(candidates)
        chosen = []
        used_sites = set()
        for interstitial_site, vacancy_site in candidates:
            if (
                interstitial_site in used_sites
                or vacancy_site in used_sites
            ):
                continue
            chosen.append((interstitial_site, vacancy_site))
            used_sites.update((interstitial_site, vacancy_site))
            if len(chosen) == num:
                break
        if len(chosen) < num:
            raise ValueError(
                f"Requested {num} Frenkel pairs, but only {len(chosen)} "
                "disjoint first-neighbor pairs are available."
            )
        selected, ghost = zip(*chosen)
        return np.asarray(selected, dtype=int), np.asarray(ghost, dtype=int)

    def action(self,system,**kwargs):
        num = kwargs['num']
        chosen, chosen_viz = self.select(system, num=num)
        for selected, ghost in zip(chosen, chosen_viz):
            self.create(system, selected, ghost)
