import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import animation
import os
from kmc_classes import *
import sys
import warnings
import PARAM

warnings.filterwarnings("ignore")   
plt.rcParams['animation.ffmpeg_path'] = r'C:\ffmpeg\ffmpeg.exe'  
identifier = 'New'

animation_mode = PARAM.animation_mode
time_limit = np.inf

#getting parameters from
rounds        = PARAM.rounds
processes     = PARAM.processes
monomolecular = PARAM.monomolecular
s1s           = PARAM.s1s
t1s           = PARAM.t1s
num_ex        = PARAM.num_ex
anni          = PARAM.anni

#reading the lattice
X	      = PARAM.X
Y             = PARAM.Y
Z             = PARAM.Z
Mats          = PARAM.Mats

#functions from morphology ( generate particles and distribuition of energies)
gen_particles   = PARAM.gen_function
energy_fun      = PARAM.ener_function

# the respective parameters of the functions above
param_genfunc  = PARAM.parameters_genfunc
param_energy   = PARAM.parameters_enefunc

anni_funcs_array = PARAM.annihi_funcs_array

           
# runs the annihilations defined in anni_funcs_array                 
def anni_general(system,Ss,anni_funcs_array):  
    mapa_singlet = []
    mapa = []
    locs = np.array([s.location() for s in Ss])
    
    if len(locs) > len(set(locs)):
        locs2 = np.array(list(set(locs)))
        for i in range(len(locs2)):
            indices = np.where(locs == locs2[i])
            if len(indices[0]) > 1:
                tipos = [Ss[j].kind() for j in indices[0]]
                
                #running all the choosen annifuncs from morphology.py
                for anni_func in anni_funcs_array:
                    anni_func(system,tipos,Ss,indices)
              	

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
            #anni_singlet(system,Ss)
            anni_general(system,Ss,anni_funcs_array)
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
        
def animate(num,system,ax): 
    Ss = step(system)
    X, Y, Z = system.get_XYZ()
    mats = system.get_mats()
    nx,ny = [],[]
    #plt.cla()
    ax.clear()
    
    X0 = X[mats == 0]
    Y0 = Y[mats == 0]
    Z0 = Z[mats == 0]
    
    X1 = X[mats == 1]
    Y1 = Y[mats == 1]
    Z1 = Z[mats == 1]
    
    ax.scatter(X0,Y0,Z0,alpha=0.2,color='black')
    ax.scatter(X1,Y1,Z1,alpha=0.2,color='blue')
    try:  
        for s in Ss:
            xs = X[s.location()]        	
            ys = Y[s.location()]
            zs = Z[s.location()]
              
            #opcao 1 (com os markers)
            '''        
            if s.kind() == 'electron':
                ax.scatter(xs,ys,zs,marker="$e^-$",color='red',s=200,alpha=1,label=s.kind())
            elif s.kind() == 'hole':
                ax.scatter(xs,ys,zs,marker="$h^+$",color='blue',s=200,alpha=1,label=s.kind())
            elif s.kind() == 'triplet':
                ax.scatter(xs,ys,zs,marker='$T_1$',color='green',s=200,alpha=1,label=s.kind())
            elif s.kind() == 'singlet':
                ax.scatter(xs,ys,zs,marker='$S_1$',color='orange',s=200,alpha=1,label=s.kind())
            ''' 
            #opcao 2(sem marker)                
            if s.kind() == 'electron':
                ax.scatter(xs,ys,zs,color='red',s=100,alpha=1,label=s.kind())
            elif s.kind() == 'hole':
                ax.scatter(xs,ys,zs,color='blue',s=100,alpha=1,label=s.kind())
            elif s.kind() == 'triplet':
                ax.scatter(xs,ys,zs,color='green',s=100,alpha=1,label=s.kind())
            elif s.kind() == 'singlet':
                ax.scatter(xs,ys,zs,color='orange',s=100,alpha=1,label=s.kind())
    except:
        pass
    
    #remove duplicates on legend    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    
    return ax,



s1, t1 = energy_fun(param_energy)


if animation_mode:
    system = System(X,Y,Z,Mats)
    system.set_s1(s1)
    system.set_t1(t1)
    system.set_orbital(t1,s1)
    excitons = gen_particles(param_genfunc)
    system.set_particles(excitons)
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    ani = animation.FuncAnimation(fig, animate, fargs=[system,ax],
                                interval=200, blit=False,repeat=False,cache_frame_data=True)#,save_count=1000) 
    #ani.save('charges.avi', fps=20, dpi=300)
    #os.system("C:\ffmpeg\ffmpeg.exe -i charges.avi charges.gif")
    
    plt.show()
else:
    for i in range(rounds):
        system = System(X,Y,Z,Mats)
        system.set_s1(s1)
        system.set_t1(t1)
        system.set_orbital(t1,s1)
        excitons = gen_particles(param_genfunc)
        system.set_particles(excitons)
        step(system)
        spectra(system)
    
            
            
        
    
