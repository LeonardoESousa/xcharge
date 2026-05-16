
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
        
        DEe = lumos - (homos[local] + s1s[local])
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
        
        DEh = (lumos[local] - s1s[local]) - homos  
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


def corrected_energies(system,s,r,dx,dy,dz):
    potential = system.electrostatic() - (
        system.field[0] * dx + system.field[1] * dy + system.field[2] * dz
    )
    r_safe = np.where(r != 0, r, np.inf)
    potential -= s.charge * abs(e) / (4 * np.pi * system.epsilon * r_safe)
    indices_e = np.fromiter(
        (pos for pos in system.positions_by_species.get("electron", set()) if pos != s.position),
        dtype=int,
    )
    indices_h = np.fromiter(
        (pos for pos in system.positions_by_species.get("hole", set()) if pos != s.position),
        dtype=int,
    )

    if s.species == 'electron':
        engs  = np.array(system.lumo, copy=True)
        if indices_h.size:
            engs[indices_h] = system.homo[indices_h]
        if indices_e.size:
            engs[indices_e] = np.inf
        engs += -1*potential
        DE = (engs - engs[s.position]) + abs(engs - engs[s.position])
    elif s.species == 'hole':
        engs  = np.array(system.homo, copy=True)
        if indices_e.size:
            engs[indices_e] = system.lumo[indices_e]
        if indices_h.size:
            engs[indices_h] = -np.inf
        engs += -1*potential
        DE = (engs[s.position] - engs) + abs(engs[s.position] - engs)
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
        mats   = kwargs['mats']
        mat    = kwargs['matlocal']
        num      = len(mats)        
        
        AtH        = raios(len(r),self.AtH,mat,self.inv,mats)
        in_loc_rad = self.inv[mat]

        DE = corrected_energies(system,particle,r,dx,dy,dz)                            	               
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
  # Initialize the Raios array with the value of Rf[(mat,mat)]
  taxas = np.empty(num)
  ks[(999,999)] =0
  taxas.fill(ks[(mat,mat)])
  
  # Use NumPy's where function to set the values of Raios for the other materials
  taxas = np.where(mats == 999, 0, taxas) # negating hop to ghost site
  for m in materials_list:
    ks[(999,m)] = 0
    if m != mat:
      taxas = np.where(mats == m, ks[(mat,m)], taxas) # where mats == m, rate will be ks[(m,m)]). To be stored in the already initialized taxas vec
  return taxas


##DEFECTS MIGRATION RATE##############################################################
class Migration:
    def __init__(self,**kwargs):
        self.kind = 'jump'
        self.k = kwargs['k']
        self.materials_list = list(set([key[0] for key in self.k.keys()]))

    def rate(self,**kwargs):
        r        = kwargs['r']
        system   = kwargs['system']
        particle = kwargs['particle']
        mats     = kwargs['mats']    
        local    = np.argwhere(r == 0)[0][0]
        mat      = mats[local]
        cut      = kwargs['cut']   # full-lattice indices of those neighbors
        
        if cut is None:
            # Fallback: assume contiguous slice from 0
            cut = np.arange(len(r))        
        
        #migration can not occur in occupied sites!
        forbidden_sites  = set().union(*(system.positions_by_species.values()))
        occupied = np.fromiter((site in forbidden_sites for site in cut),
                               dtype=bool, count=len(r))
                               
                               
        taxa = filter(len(mats),self.k,mat,mats,self.materials_list)
        taxa = np.where(occupied, 0.0, taxa) # if occupied, rate is 0, else, remain the calc value

        taxa[local] = 0
        
        return taxa


    def action(self,particle,system,local):
        particle.move(local,system)
class DissociationFP:
    def __init__(self,**kwargs):
        self.kind = 'dissociation_fp'
        self.k = kwargs['k']
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
        self.k = kwargs['k']
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
        self.k = kwargs['k']
        self.materials_list = list(set([key[0] for key in self.k.keys()]))
    def rate(self, **kwargs):
        r        = kwargs['r']
        system   = kwargs['system']
        particle = kwargs['particle']
        mats     = kwargs['mats']    
        local    = np.argwhere(r == 0)[0][0]
        mat      = mats[local]
        
        cut      = kwargs['cut']   # full-lattice indices of those neighbors
        
        if cut is None:
            # Fallback: assume contiguous slice from 0
            cut = np.arange(len(r))
            
        sites_above_FPradius  = set(np.where(r > system.FPcutoff)[0])  #formation/generation cutoff should be more severe than hopping cutoff
        
        origin_type = particle.species
        #two forms of formation: I first  or V first
        if origin_type == 'interstitial':
            neigh_type = 'vacancy'
        if origin_type == 'vacancy':
            neigh_type = 'interstitial'        
        
        interstitial_sites = system.positions_by_species.get(neigh_type, set())
        interstitial_sites = interstitial_sites- sites_above_FPradius
        # Boolean occupancy mask aligned with r/mats via cut
        occupied = np.fromiter((site in interstitial_sites for site in cut),
                               dtype=bool, count=len(r))

        # Rate: k when destination has an insterstitial, else 0
        taxa = filter(len(mats),self.k,mat,mats,self.materials_list)
        taxa = np.where(occupied, taxa, 0.0)

        # Never hop to self
        taxa[r == 0] = 0
        return taxa


    def action(self, particle, system, local):
        # FP generation requires the addition of a ghost(deactivated) site
        # where new particles cant hop to it 
        # FP.ghost_site is the position (in len(system.X))
        # material type 999 is reserved for this special site
        #convention:
        # ->local where the FP is centered, conventioned to be the place where the Intersitial moves to
        # ->ghost position is reserved for the vacancy
        FP = Frenkelpair(local)
        FP.ghost_site = particle.position
        FP.origin_site = system.mats[particle.position]
        system.mats[particle.position] = 999
        system.set_particles([FP])
        particle.kill(self.kind, system, system.s1, 'converted')
        #debug
        #xfp,yfp,zfp = system.X[particle.position],system.Y[particle.position],system.Z[particle.position]
        #xg,yg,zg = system.X[local],system.Y[local],system.Z[local]
        #dx,dy,dz =  xg-xfp, yg-yfp, zg-zfp
        #dist = np.sqrt(dx**2 + dy**2+dz**2)
        #print('form:',dist,particle.position,local)
        inters = [ p for p in system.particles if p.species=="interstitial"]
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
'''
class FP_generation:
    def __init__(self,**kwargs):
        self.kind = 'generation'
        self.k = kwargs['k']
        
    def rate(self,**kwargs):
        return self.k#[kwargs['material']]
     
    def action(self,system,**kwargs):
        npart            = kwargs['npart']
        pairs            = kwargs['pairs']        
        
        forbidden_sites  = set().union(*(system.positions_by_species.values()))
        all_sites        = set(range(len(system.X)))
        sites999         = np.where(system.mats == 999)[0] #cant create FP at 999 materials
        avail_sites      = list(all_sites-forbidden_sites-set(sites999))
        avail_sites      = random.sample(avail_sites, k=len(avail_sites))        
        for i, site in enumerate(avail_sites):
            orig_mat     = system.mats[site]
            dx, dy, dz   = system.distance(system,site)
            hopsites     = avail_sites.copy()
            r            = kmc.utils.distances(dx,dy,dz,len(dx))
            uncut        = np.where(r > system.FPcutoff)[0] 
            mat_mismatch = np.where(system.mats == orig_mat)[0] #FP forms occupying two sites of different composition
            hopsites     = list(set(hopsites) - set(uncut) - set(mat_mismatch) - set([i]) -set(sites999))
            if len(hopsites) > 0:
                selec = i
                viz   = hopsites
                break
        selected = avail_sites[selec]
        ghost    = random.sample(viz, 1)[0]
        xfp,yfp,zfp = system.X[selected],system.Y[selected],system.Z[selected]
        xg,yg,zg = system.X[ghost],system.Y[ghost],system.Z[ghost]
        dx,dy,dz =  xg-xfp, yg-yfp, zg-zfp
        dist = np.sqrt(dx**2 + dy**2+dz**2)
        print('gen:',dist)
        FP = Frenkelpair(selected)
        FP.ghost_site = ghost
        FP.origin_site = system.mats[ghost] #storing the old material type
        system.mats[ghost] = 999 #changing the mat
        system.set_particles([FP])
        FP.make_text(system,system.s1,causamortis='generated')
        FP.stamp_time(system)
'''
####################################################################################
class FP_generation:
    def __init__(self,**kwargs):
        self.kind = 'generation'
        self.k = kwargs['k']
        self.pairs = kwargs['pairs']
        try:
          self.causamortis = kwargs['causamortis']
        except:
          self.causamortis = "generated"
    def rate(self,**kwargs):
        return self.k#[kwargs['material']]
    def create(self,system,selected,ghost,**kwargs):
        FP = Frenkelpair(selected)
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
        pairs            = self.pairs
        chosen,chosen_viz = [],[]
        for pair in pairs:
            V_site_type,I_site_type = pair        
            forbidden_sites  = set().union(*(system.positions_by_species.values()))
            all_I_sites      = np.where(system.mats == I_site_type)[0]
            all_V_sites      = np.where(system.mats == V_site_type)[0]
            sites999         = np.where(system.mats == 999)[0] #cant create FP at 999 materials
            avail_sites      = list(set(all_I_sites)-forbidden_sites-set(sites999))
            avail_sites      = random.sample(avail_sites, k=len(avail_sites))        
            for i, site in enumerate(avail_sites):
              orig_mat     = system.mats[site]
              dx, dy, dz   = system.distance(system,site)
              r            = kmc.utils.distances(dx,dy,dz,len(dx))
              uncut        = np.where(r > system.FPcutoff)[0] 
              hopsites     = list(set(all_V_sites) - set(uncut) - set([i]) -set(sites999) -set(forbidden_sites))
              if len(hopsites) > 0:
                  selec = i
                  viz   = hopsites
                  break
            chosen.append(avail_sites[selec])
            chosen_viz.append(random.sample(viz, 1)[0])
        return chosen,chosen_viz        
    def action(self,system,**kwargs):
        try:
          chosen,chosen_viz = self.select(system)
          for i, (selec,viz) in enumerate(zip(chosen,chosen_viz)):
            self.create(system,selec,viz)
        except Exception as e:
          print(f'This is a warning, generation/creation was not successful for some reason! {e}')

