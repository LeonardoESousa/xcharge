import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import animation
import os
from kmc_classes import *
import sys
import warnings
import PARAM
import morphology

warnings.filterwarnings("ignore")   
plt.rcParams['animation.ffmpeg_path'] = r'C:\ffmpeg\ffmpeg.exe'  
identifier = 'New'

animation_mode = True
time_limit = np.inf


rounds        = PARAM.rounds
processes     = PARAM.processes
monomolecular = PARAM.monomolecular
s1s           = PARAM.s1s
t1s           = PARAM.t1s
num_ex        = PARAM.num_ex
anni          = PARAM.anni

gen_singlets   = morphology.generate_function
homo_lumo      = morphology.energy_function

def read_lattice(file_name):
	X,Y,Z,Mats = [], [], [], []
	
	with open(file_name, 'r') as f:
		for line in f:
			line = line.split()
			
			x    = float(line[0])
			y    = float(line[1])
			z    = float(line[2])
			mat  = int(float(line[3]))
			
			
			X.append(x)
			Y.append(y)
			Z.append(z)
			Mats.append(mat)
	X = np.array(X)
	Y = np.array(Y)	
	Z = np.array(Z)	
	Mats = np.array(Mats)
	return X,Y,Z,Mats
	
	



def anni_singlet(system,Ss):  
    mapa_singlet = []
    mapa = []
    locs = np.array([s.location() for s in Ss])
    if len(locs) > len(set(locs)):
        locs2 = np.array(list(set(locs)))
        for i in range(len(locs2)):
            indices = np.where(locs == locs2[i])
            if len(indices[0]) > 1:
                tipos = [Ss[j].kind() for j in indices[0]]
                if 'electron' in tipos and 'hole' in tipos:        
                    Ss[indices[0][tipos.index('electron')]].kill('anni',system,system.get_s1())
                    Ss[indices[0][tipos.index('hole')]].kill('anni',system,system.get_s1())
                    if random.uniform(0,1) <= 0.75:
                        system.add_particle(Exciton('triplet',locs[indices[0][0]]))
                    else:
                        system.add_particle(Exciton('singlet',locs[indices[0][0]]))

    
    #for s in Ss:
    #    if s.kind() != 'singlet':
    #        loc = s.location()
    #        if loc in mapa_singlet:
    #            s.kill('anni',system,system.get_s1())
    #        else:
    #            mapa_singlet.append(s.location())

def decision(s,system):
    kind = s.kind()
    local = s.location()
    X,Y,Z = system.get_XYZ()
    Mat   = system.get_mats()[local]
    dx = np.nan_to_num(X - X[local]) 
    dy = np.nan_to_num(Y - Y[local])
    dz = np.nan_to_num(Z - Z[local])
    r  = np.sqrt(dx**2+dy**2+dz**2)
    r[r == 0] = np.inf
    
    final_rate = []
    labels     = []
    chosen     = []
    hop = processes.get(kind)
    for transfer in hop:    
        jump_rate  = transfer.rate(r=r,system=system,particle=s)
        probs = np.cumsum(jump_rate)/np.sum(jump_rate)
        sorte = random.uniform(0,1)
        try:
            chosen.append(np.where(sorte < probs)[0][0])
        except:
            chosen.append(local)
        final_rate.append(jump_rate[chosen[-1]])
        labels.append(transfer)
 
    mono = monomolecular.get(kind)   
    for m in mono:
        final_rate.append(m.rate(material=Mat))
        labels.append(m)
        chosen.append(local)
    
    probs = np.cumsum(final_rate)/np.sum(final_rate)
    sorte = random.uniform(0,1)
    jump = np.where(sorte < probs)[0][0]
    dt = (1/np.sum(final_rate))*np.log(1/random.uniform(1E-12,1))
    #print(kind,Mats[local],Mats[chosen],probs,labels[jump],dt)
    return labels[jump], chosen[jump], dt
 
  
def step(system): 
    while system.count_particles() > 0 and system.clock() < time_limit:
        Xx = system.get_particles()
        Ss = list(Xx)
        #print([m.kind() for m in Ss])
        J, W, DT = [],[],[]
        for s in Ss:
            jump, where, dt = decision(s,system)
            J.append(jump)
            W.append(where)
            DT.append(dt)    
        #print(W)
        time_step = min(DT)
        realtime = system.clock() + time_step
        system.set_clock(realtime)
        fator = random.uniform(0,1)
        for i in range(len(Ss)):
            if fator <= time_step/DT[i]:
                J[i].action(Ss[i],system,W[i])
        if anni:
            anni_singlet(system,Ss)
        if animation_mode:
            return Ss       
    Xx = system.get_particles()
    Ss = list(Xx)
    for s in Ss:
        Ss[i].kill('alive',system,system.get_s1())
  
            
                
#Prints Spectra
def spectra(system):
    if os.path.isfile("Simulation_"+identifier+".txt") == False:
        with open("Simulation_"+identifier+".txt", "w") as f:
            texto = "{0:^10}    {1:^6} {2:^6} {3:^6} {4:^4} {5:^4} {6:^9} {7:^6} {8:^6} {9:^6} {10:^4}".format("Time", "DeltaX", "DeltaY", "DeltaZ", "Type", "Energy", "Location" ,"FinalX", "FinalY", "FinalZ", "Causa Mortis")
            f.write(texto+"\n") 
    with open("Simulation_"+identifier+".txt", "a") as f:   
        for s in system.get_dead():
            texto = s.write()
            f.write(texto)
        f.write("Fim\n")
        
def animate(num,system,line): 
    Ss = step(system)
    X, Y, Z = system.get_XYZ()
    mats = system.get_mats()
    nx,ny = [],[]
    plt.cla()
    line.axis('off')
    X0 = X[mats == 0]
    Y0 = Y[mats == 0]
    
    
    X1 = X[mats == 1]
    Y1 = Y[mats == 1]
    
    
    line.plot(X0,Y0,'s',color='black',markersize=2)
    line.plot(X1,Y1,'s',color='blue' ,markersize=2)
    try:
        for s in Ss:
            line.plot(X[s.location()],Y[s.location()],marker="s",color='white', markersize=12)
            if s.kind() == 'electron':
                line.plot(X[s.location()],Y[s.location()],marker="$e^-$",color='red', markersize=12)
            elif s.kind() == 'hole':
                line.plot(X[s.location()],Y[s.location()],marker="$h^+$",color='blue', markersize=12)
            elif s.kind() == 'triplet':
                line.plot(X[s.location()],Y[s.location()],marker='$T_1$',color='green', markersize=12)
            elif s.kind() == 'singlet':
                line.plot(X[s.location()],Y[s.location()],marker='$S_1$',color='orange', markersize=12)
            #nx.append(X[s.location()])
            #ny.append(Y[s.location()])
    except:
        pass
    return line,


X,Y,Z,Mats = read_lattice("lattice.txt")

s1, t1 = homo_lumo(s1s, t1s, Mats)


if animation_mode:
    system = System(X,Y,Z,Mats)
    system.set_s1(s1)
    system.set_t1(t1)
    system.set_orbital(t1,s1)
    excitons = gen_singlets(num_ex,len(X))
    system.set_particles(excitons)
    fig, line = plt.subplots()
    line.axis('off')
    ani = animation.FuncAnimation(fig, animate, fargs=[system,line],
                                interval=200, blit=False,repeat=False)#,save_count=1000) 
    #ani.save('charges.avi', fps=20, dpi=300)
    #os.system("C:\ffmpeg\ffmpeg.exe -i charges.avi charges.gif")
    plt.show()
else:
    for i in range(rounds):
        system = System(X,Y,Z,Mats)
        system.set_s1(s1)
        system.set_t1(t1)
        system.set_orbital(t1,s1)
        excitons = gen_singlets(num_ex,len(X))
        system.set_particles(excitons)
        step(system)
        spectra(system)
    
            
            
        
    
