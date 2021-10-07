import numpy as np
import random
from kmc.kmc_classes import *
import sys
import warnings
import os
import matplotlib.pyplot as plt
from matplotlib import animation
from joblib import Parallel, delayed
from mpl_toolkits.mplot3d import Axes3D
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

    try:
    	jump = np.where(sorte < probs)[0][0]
    	dt = (1/np.sum(final_rate))*np.log(1/random.uniform(1E-12,1))
    #If no option is available, particle stands still
    except:
    	jump = 0
    	dt = np.inf
    
    return labels[jump], chosen[jump], dt

def step(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        Ss = system.particles.copy()     
        J, W, DT = [],[],[]
        for s in Ss:
            jump, where, dt = decision(s,system)
            J.append(jump)
            W.append(where)
            DT.append(dt)    
        time_step = min(DT)
        system.time += time_step
        fator = random.uniform(0,1)
        for i in range(len(Ss)):
            if fator <= time_step/DT[i]:
                J[i].action(Ss[i],system,W[i])
        if system.anni:
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
        
def animate(num,system,ax,marker_option): 
    Ss = step(system)
    X, Y, Z = system.X, system.Y, system.Z        
    mats = system.mats                            
    ax.clear()
    #ploting the sites according to mat index
    colors_dic = {0:'black', 1:'blue', 2:'red', 3:'green', 4:'yellow'}
    n_mats = np.unique(mats)
    for mat in n_mats:
        X_mat = X[mats == mat]
        Y_mat = Y[mats == mat]
        Z_mat = Z[mats == mat]
        ax.scatter(X_mat,Y_mat,Z_mat,alpha=0.25,color=colors_dic.get(int(mat)))
    try:  
        for s in Ss:
            xs = X[s.position]        	
            ys = Y[s.position]
            zs = Z[s.position]    
                
            if marker_option == 1:      
                if s.species == 'electron':
                    ax.scatter(xs,ys,zs,marker="$e^-$",color='red',s=200,alpha=1,label=s.species)
                elif s.species == 'hole':
                    ax.scatter(xs,ys,zs,marker="$h^+$",color='blue',s=200,alpha=1,label=s.species)
                elif s.species == 'triplet':
                    ax.scatter(xs,ys,zs,marker='$T_1$',color='green',s=200,alpha=1,label=s.species)
                elif s.species == 'singlet':
                    ax.scatter(xs,ys,zs,marker='$S_1$',color='orange',s=200,alpha=1,label=s.species)
            if marker_option == 0:               
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
    ax.text2D(0.03, 0.94, "eps  = %.2f"    % (system.eps_rel), transform=ax.transAxes) #eps
    ax.text2D(0.03, 0.90, "npart  = %.0f"  % (len(system.particles)), transform=ax.transAxes) #npart
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')     

    #pausing in the first frame
    #if system.pause:
    #    ani.event_source.stop()
    return ax,

def draw_lattice(X,Y,Z,Mats,color_dir,fig_name):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(X, Y, Z,c=color_dir,marker='^');
    
    try:
        plt.show()    
    except:
        plt.switch_backend('agg')
        plt.savefig(fig_name+'.png')
        
#RUN of a single round   
def RUN(system):
    step(system)
    spectra(system)
#n_rounds_runned = 0   
def make_system(module_param):

    #getting the lattice for this particular round (can be a fresh one or not, depending on the user's choice)
    X,Y,Z,Mats = module_param.lattice_func(module_param.lattice_func_par)
    system = System(X,Y,Z,Mats)   
    
    system.set_basic_info(module_param.monomolecular,module_param.processes,
    module_param.identifier,module_param.animation_mode,module_param.time_limit,module_param.pause,
    module_param.anni,module_param.annihi_funcs_array) 
    
    #getting the s1, t1, HOMO, and LUMO and then adding to the system
    s1,t1,HOMO,LUMO = module_param.ener_function(module_param.parameters_enefunc,Mats)
    system.set_energies(s1,t1,HOMO,LUMO)
    
    try:
        system.set_dipoles(module_param.dipoles)
    except:
        pass
    system.set_medium(module_param.relative_eps)
    
    #setting up particle generation
    selection          = range(len(X))
    parameters_genfunc = [module_param.num_ex,selection]
    
    excitons = param.gen_function(parameters_genfunc)
    system.set_particles(excitons)
    
    return system 
    
#setting up the animation object and adding responses to events    
def run_animation(param):
    ani_running = True

    def onClick(event): #if somenone clicks on the ani, this happens
        nonlocal ani_running
        if ani_running:
            ani.event_source.stop()
            ani_running = False
        else:
            ani.event_source.start()
            ani_running = True

    def pause_plot(event,pause): #if pause = true, this will happen
        nonlocal ani_running
        if pause:
            ani.event_source.stop()
            ani_running  = False

    system = make_system(param)
    pause           = param.pause # to freeze on the first frame
    marker_type     = param.marker_type
                    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fig.canvas.mpl_connect('button_press_event', onClick) #pausing if clicking
    fig.canvas.mpl_connect('draw_event', lambda event: pause_plot(event, pause)) #pausing if pause = True at the first frame
   
    ani = animation.FuncAnimation(fig, animate, fargs=[system,ax,marker_type],
                                    interval=25, blit=False,repeat=False,cache_frame_data=True)#,save_count=1000)  
                             
                                       
    return ani 
    

    
def main():

    n_proc          = param.n_proc
    identifier      = param.identifier 
    animation_mode  = param.animation_mode
    save_animation  = param.save_animation 
    animation_exten = param.animation_exten    
    rounds          = param.rounds
    
    if animation_mode:
        ani = run_animation(param)
        path=identifier+"_animation."+animation_exten
                                                   
        if save_animation:                   
            
            #save .gif
            if animation_exten == 'gif':
                ani.save(path, writer='imagemagick', fps=10)
            
            #save .mp4
            if animation_exten == 'mp4':
                writervideo = animation.FFMpegWriter(fps=10) 
                ani.save(path, writer=writervideo)
        
        plt.show()
                
    else:
    
        Parallel(n_jobs=n_proc, backend = 'loky')(delayed(RUN)(make_system(param)) for _ in range(rounds))
                
if __name__ == "__main__":
    sys.exit(main())        
