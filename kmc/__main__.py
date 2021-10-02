import numpy as np
import random
from kmc.kmc_classes import *
import sys
import warnings
import os
import matplotlib.pyplot as plt
from matplotlib import animation
from joblib import Parallel, delayed
import importlib
warnings.filterwarnings("ignore")   


working_dir = os.getcwd()+'/'
#importing param module

spec  = importlib.util.spec_from_file_location(sys.argv[1].split('.')[0], working_dir+sys.argv[1])
param = importlib.util.module_from_spec(spec)
spec.loader.exec_module(param)

         

# runs the annihilations defined in anni_funcs_array                 
def anni_general(system,Ss,anni_funcs_array):   
    locs = np.array([s.position for s in Ss])
    
    if len(locs) > len(set(locs)):
        locs2 = np.array(list(set(locs)))
        for i in range(len(locs2)):
            indices = np.where(locs == locs2[i])
            if len(indices[0]) > 1:

                tipos = [Ss[j].species for j in indices[0]]
                
                #running all the choosen annifuncs from morphology.py
                for anni_func in anni_funcs_array:
                    anni_func(system,tipos,Ss,indices,locs)
              	


def decision(s,system):
    kind = s.species      
    local = s.position    
    X,Y,Z = system.X, system.Y, system.Z 
    Mat   = system.mats[local]   
    dx = np.nan_to_num(X - X[local]) 
    dy = np.nan_to_num(Y - Y[local])
    dz = np.nan_to_num(Z - Z[local])
    r  = np.sqrt(dx**2+dy**2+dz**2)
    r[r == 0] = np.inf
    
    final_rate = []
    labels     = []
    chosen     = []
    hop = system.processes.get(kind)
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
 
    mono = system.monomolecular.get(kind)   
    for m in mono:
        final_rate.append(m.rate(material=Mat))
        labels.append(m)
        chosen.append(local)
    
    probs = np.cumsum(final_rate)/np.sum(final_rate)
    sorte = random.uniform(0,1)
    
    '''
    jump = np.where(sorte < probs)[0][0]
    dt = (1/np.sum(final_rate))*np.log(1/random.uniform(1E-12,1))
    '''
    
    try:
    	jump = np.where(sorte < probs)[0][0]
    	dt = (1/np.sum(final_rate))*np.log(1/random.uniform(1E-12,1))
    except:
    	jump = 0
    	dt = np.inf
    
    #print(dt)
    #input()	
    #print(kind,Mats[local],Mats[chosen],probs,labels[jump],dt)
    return labels[jump], chosen[jump], dt
 
#checks if a electron's jumb matches with a hole position or vice versa 
def pair_matching(particles_sample,W):
    parts = particles_sample.copy()
    W.copy()
    pos   = [part.position for part in parts ]
    jumps = W.copy()
    
    n = len(parts)
    
    recomb = []
    
    for i in range(n):
    	part_i = parts[i]
    	if part_i.species == 'hole':
    	    for j in range(n):
                part_j = parts[j]
                if part_j.species == 'electron':
                    
                    #print(pos[i],jumps[i],"    ",pos[j],jumps[j])
                    
                    if jumps[i] == pos[j] or jumps[j] == pos[i]: #or jumps[i] == jumps[j]:
                        recomb.append( i )
                        
    return recomb
            
 

def step(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        Ss = system.particles.copy()     
        #print([m.species for m in Ss])
        J, W, DT = [],[],[]
        for s in Ss:
            jump, where, dt = decision(s,system)
            J.append(jump)
            W.append(where)
            DT.append(dt)    
        #print(W)
        time_step = min(DT)
        system.time += time_step
        #print(system.time) 
        fator = random.uniform(0,1)
        for i in range(len(Ss)):
            if fator <= time_step/DT[i]:
                J[i].action(Ss[i],system,W[i])
        if system.anni:
            #anni_singlet(system,Ss)
            anni_general(system,Ss,system.anni_funcs_array)
        if system.animation_mode:
            return Ss       
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1)
  
  
            
                
#Prints Spectra
def spectra(system):
    if os.path.isfile("Simulation_"+system.identifier+".txt") == False:
        with open("Simulation_"+system.identifier+".txt", "w") as f:
            texto = "{0:^10}    {1:^6} {2:^6} {3:^6} {4:^4} {5:^4} {6:^9} {7:^6} {8:^6} {9:^6} {10:^4}".format("Time", "DeltaX", "DeltaY", "DeltaZ", "Type", "Energy", "Location" ,"FinalX", "FinalY", "FinalZ", "Causa Mortis")
            f.write(texto+"\n") 
    with open("Simulation_"+system.identifier+".txt", "a") as f:   
        for s in system.dead:
            texto = s.write()
            f.write(texto)
        f.write("Fim\n")
        
def animate(num,system,ax): 
    Ss = step(system)
    X, Y, Z = system.X, system.Y, system.Z        
    mats = system.mats                            
    ax.clear()
    
    X0 = X[mats == 0]
    Y0 = Y[mats == 0]
    Z0 = Z[mats == 0]
    
    X1 = X[mats == 1]
    Y1 = Y[mats == 1]
    Z1 = Z[mats == 1]
    
    #printing the lattice
    ax.scatter(X0,Y0,Z0,alpha=0.1,color='black')
    ax.scatter(X1,Y1,Z1,alpha=0.1,color='blue')
    
    try:  
        for s in Ss:
            xs = X[s.position]        	
            ys = Y[s.position]
            zs = Z[s.position]
              
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
            if s.species == 'electron':
                ax.scatter(xs,ys,zs,color='red',s=100,alpha=1,label=s.species)
            elif s.species == 'hole':
                ax.scatter(xs,ys,zs,color='blue',s=100,alpha=1,label=s.species)
            elif s.species == 'triplet':
                ax.scatter(xs,ys,zs,color='green',s=100,alpha=1,label=s.species)
            elif s.species == 'singlet':
                ax.scatter(xs,ys,zs,color='orange',s=100,alpha=1,label=s.species)
    except:
        pass
    
    #remove duplicates on legend    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    ax.text2D(0.03, 0.98, "time = %.2e ps" % (system.time), transform=ax.transAxes) #time
    ax.text2D(0.03, 0.94, "eps  = %.2f" %    (system.eps_rel), transform=ax.transAxes) #eps
    ax.text2D(0.03, 0.90, "npart  = %.0f" % (len(system.particles)), transform=ax.transAxes) #npart
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')     

    #pausing in the first frame
    if system.pause:
        ani.event_source.stop()
    return ax,


#RUN of a single round   
def RUN(system):
    step(system)
    spectra(system)
   

def main():
    n_proc          = param.n_proc
    identifier      = param.identifier 
    animation_mode  = param.animation_mode
    time_limit      = param.time_limit 
    pause           = param.pause # to freeze on the first frame

    #getting parameters from
    rounds           = param.rounds
    processes        = param.processes
    monomolecular    = param.monomolecular
    anni             = param.anni
    anni_funcs_array = param.annihi_funcs_array

    if animation_mode:
    
        #path="animation.mp4"
            
    
        system = param.make_system()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
           
    
        ani = animation.FuncAnimation(fig, animate, fargs=[system,ax],
                                    interval=25, blit=False,repeat=False,cache_frame_data=True)#,save_count=1000) 
        #ani.save('charges.avi', fps=20, dpi=300)
        #os.system("C:\ffmpeg\ffmpeg.exe -i charges.avi charges.gif")
        
        #salvar .gif
        #ani.save(path, writer='imagemagick', fps=30)
        
        #salvar .mp4
        #writervideo = animation.FFMpegWriter(fps=10) 
        #ani.save(path, writer=writervideo)
        
        plt.show()
        
        
    else:
    
        Parallel(n_jobs=n_proc, backend = 'loky')(delayed(RUN)(param.make_system()) for i in range(rounds))
          

if __name__ == "__main__":
    sys.exit(main())        
