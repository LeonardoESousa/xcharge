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

#usuario de ruindows
#plt.rcParams['animation.ffmpeg_path'] = r'C:\ffmpeg\ffmpeg.exe'  

identifier     = PARAM.identifier 
animation_mode = PARAM.animation_mode
time_limit     = PARAM.time_limit 
pause          = PARAM.pause # to freeze on the first frame

#getting parameters from
rounds        = PARAM.rounds
processes     = PARAM.processes
monomolecular = PARAM.monomolecular
anni          = PARAM.anni

anni_funcs_array = PARAM.annihi_funcs_array

  

         
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
 
  
def step(system): 
    while system.count_particles() > 0 and system.time < time_limit:
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
        if anni:
            #anni_singlet(system,Ss)
            anni_general(system,Ss,anni_funcs_array)
        if animation_mode:
            return Ss       
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1)
  
            
                
#Prints Spectra
def spectra(system):
    if os.path.isfile("Simulation_"+identifier+".txt") == False:
        with open("Simulation_"+identifier+".txt", "w") as f:
            texto = "{0:^10}    {1:^6} {2:^6} {3:^6} {4:^4} {5:^4} {6:^9} {7:^6} {8:^6} {9:^6} {10:^4}".format("Time", "DeltaX", "DeltaY", "DeltaZ", "Type", "Energy", "Location" ,"FinalX", "FinalY", "FinalZ", "Causa Mortis")
            f.write(texto+"\n") 
    with open("Simulation_"+identifier+".txt", "a") as f:   
        for s in system.dead:
            texto = s.write()
            f.write(texto)
        f.write("Fim\n")
        
def animate(num,system,ax): 
    Ss = step(system)
    X, Y, Z = system.X, system.Y, system.Z        
    mats = system.mats                            
    nx,ny = [],[]
    #plt.cla()
    ax.clear()
    

    X0 = X[mats == 0]
    Y0 = Y[mats == 0]
    Z0 = Z[mats == 0]
    
    X1 = X[mats == 1]
    Y1 = Y[mats == 1]
    Z1 = Z[mats == 1]
    
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
    return ax,


      
if animation_mode:

    #path="/home/tiago/Documents/Pesquisa/Estrutura_eletronica/KMC_TRY/KMC/animation.gif"
    path="/home/tiago/Documents/Pesquisa/Estrutura_eletronica/KMC_TRY/KMC/animation.mp4"
        

    system = PARAM.make_system()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
       
   
    ani = animation.FuncAnimation(fig, animate, fargs=[system,ax],
                                interval=25, blit=False,repeat=False,cache_frame_data=True)#,save_count=1000) 
    #ani.save('charges.avi', fps=20, dpi=300)
    #os.system("C:\ffmpeg\ffmpeg.exe -i charges.avi charges.gif")
    
    #salvar .gif
    #ani.save(path, writer='imagemagick', fps=30)
    
    #salvar .mp4
    #writervideo = animation.FFMpegWriter(fps=30) 
    #ani.save(path, writer=writervideo)
    
    plt.show()
    
else:
    for i in range(rounds):
    
        system = PARAM.make_system()
        step(system)
        spectra(system)
