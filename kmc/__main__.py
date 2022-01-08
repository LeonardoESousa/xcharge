import numpy as np
import random
from kmc.rates import *
from kmc.particles import *
from kmc.system import System
import sys
import warnings
import os
import copy
import matplotlib.pyplot as plt
from matplotlib import animation
#from joblib import Parallel, delayed
from mpl_toolkits.mplot3d import Axes3D
import importlib
warnings.filterwarnings("ignore")   
from tqdm.contrib.concurrent import thread_map, process_map


#importing param module
working_dir = os.getcwd()+'/'
spec  = importlib.util.spec_from_file_location(sys.argv[1].split('.')[0], working_dir+sys.argv[1])
param = importlib.util.module_from_spec(spec)
spec.loader.exec_module(param)

argumentos = []  
for name, value in vars(param).items():
    if hasattr(value, 'assign_to_system'):        
        argumentos.append(value)
        #print(name,hasattr(value, 'make'))

#getting all essential info from user's input
n_proc          = param.n_proc
identifier      = param.identifier 
animation_mode  = param.animation_mode
save_animation  = param.save_animation 
animation_exten = param.animation_exten    
rounds          = param.rounds
time_limit      = param.time_limit  
pause           = param.pause
marker_type     = param.marker_type
rotate          = param.rotate
frozen_lattice  = param.frozen

monomolecular       = param.monomolecular
processes           = param.processes
bimolec             = param.bimolec
bimolec_funcs_array = param.bimolec_funcs_array

def passar(*args):
    pass

def anni(system,Ss,array):
    anni_general(system,Ss,array)
   

if bimolec:
    bi_func = anni
else:
    bi_func = passar

def make_system():
    #Create instance of system
    system = System()
    #Sets system properties  
    for argumento in argumentos:
        argumento.assign_to_system(system)
    system.set_basic_info(monomolecular,processes,identifier,animation_mode,time_limit,pause,bimolec,bimolec_funcs_array) 
 
    return system 
syst = make_system() #global system-type object to be used in frozen simulations, must be kept here!
    
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
    dx = X - X[local]   
    dy = Y - Y[local]  
    dz = Z - Z[local]  
    r  = np.sqrt(dx*dx+dy*dy+dz*dz)
    
    hop = system.processes.get(kind) 
    mono = system.monomolecular.get(kind)     
    jump_rate = [transfer.rate(r=r,system=system,particle=s) for transfer in hop]
    
    try:
        #locais    = np.array([np.random.choice(np.arange(len(x)),p=x/np.sum(x)) for x in jump_rate]).astype(int)
        interval  = np.arange(len(X))
        locais    = np.array([random.choices(interval,weights=x)[0] for x in jump_rate]).astype(int)
        jump_rate = np.array([jump_rate[i][locais[i]] for i in np.arange(len(locais))])
    except:
        locais    = np.array([local])
        jump_rate = np.array([0])

    mono_rate = np.array([m.rate(material=Mat) for m in mono])
    jump_rate = np.append(jump_rate,mono_rate)
    locais2   = np.zeros(len(mono_rate)) + local
    locais    = np.append(locais,locais2.astype(int))
    labels = hop+mono 

    #jump = np.random.choice(np.arange(len(jump_rate)),p=jump_rate/np.sum(jump_rate))
    jump = random.choices(np.arange(len(jump_rate)),weights=jump_rate)[0]
    labels[jump].action(s,system,locais[jump])

    return np.sum(jump_rate)

########ITERATION FUNCTIONS#######################################################
def step_ani(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()     
        R = [decision(s,system) for s in Ss]
        system.time += (1/np.sum(R))*np.log(1/random.uniform(1E-12,1))    
        bi_func(system,Ss,system.bimolec_funcs_array)
        return Ss       
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1)
  
def step_nonani(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()     
        R = [decision(s,system) for s in Ss]
        system.time += (1/np.sum(R))*np.log(1/random.uniform(1E-12,1))    
        bi_func(system,Ss,system.bimolec_funcs_array)       
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1)
##########################################################################################

if animation_mode:
    step = step_ani
else:
    step = step_nonani


#Prints Spectra
def spectra(system):
    if os.path.isfile("Simulation_"+system.identifier+".txt") == False:
        with open("Simulation_"+system.identifier+".txt", "w") as f:
            texto = "{0:<15}  {1:<10}  {2:<10}  {3:<10}  {4:<10}  {5:<10}  {6:<10}  {7:<10}  {8:<10}  {9:<10}  {10:<10}".format("Time", "DeltaX", "DeltaY", "DeltaZ", "Type", "Energy", "Location" ,"FinalX", "FinalY", "FinalZ", "CausaMortis")
            f.write(texto+"\n") 
    with open("Simulation_"+system.identifier+".txt", "a") as f:   
        for s in system.dead:
            texto = s.write()
            f.write(texto)
        f.write("END\n")
        
def animate(num,system,ax,marker_option,rotate): 
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
                ax.scatter(xs,ys,zs,marker=s.marker,color=s.color,s=200,alpha=1,label=s.species)      
            if marker_option == 0:
                ax.scatter(xs,ys,zs,color=s.color,s=100,alpha=1,label=s.species)               
    except:
        pass
    
    if rotate:#rotating the animation by an angle of IT
        ax.view_init(azim = system.IT)
    
    #removing duplicates on legend    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    ax.text2D(0.03, 0.98, "time = %.2e ps" % (system.time), transform=ax.transAxes) #time
    ax.text2D(0.03, 0.94, "npart  = %.0f"  % (len(system.particles)), transform=ax.transAxes) #npart
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')     

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
        
#resets particles' initial position for a given system
def reroll_system(system):
    system.reset_particles()
    for argumento in argumentos:
        class_name = argumento.__class__.__name__
        if (class_name == "Create_Particles"):
            argumento.assign_to_system(system)
    
    '''
    #debug
    p  = system.particles
    pp = [ part.position for part in p ]
    print(pp,system.s1)    
    '''
    return system
            
#RUN of a single round       
def RUN(dynamic): #ROUND DYNAMICS WHERE, FOR EACH INSTANCE, THE LATTICE IS RECALCULATED
    system = make_system()
    step(system)
    spectra(system)
 
def RUN_FREEZE(dynamic): #ROUND DYNAMICS WHERE, FOR EACH INSTANCE, THE LATTICE REMAINS INTACT
    system = reroll_system(copy.deepcopy(syst))  
    step(system)
    spectra(system)

#setting up the animation object and adding responses to events    
def run_animation():
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

    system = make_system()
                    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fig.canvas.mpl_connect('button_press_event', onClick) #pausing if clicking
    fig.canvas.mpl_connect('draw_event', lambda event: pause_plot(event, pause)) #pausing if pause = True at the first frame
   
    ani = animation.FuncAnimation(fig, animate, fargs=[system,ax,marker_type,rotate],
                                    interval=25, blit=False,repeat=False,cache_frame_data=True)#,save_count=1000)  
                             
                                       
    return ani 
    

    
def main():
        
    if animation_mode:
        ani = run_animation()
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
        if not frozen_lattice: # at every round, the entire lattice is recalculated
            #Parallel(n_jobs=n_proc, backend = 'loky')(delayed(RUN)(_) for _ in range(rounds))
            process_map(RUN,range(rounds),max_workers = n_proc)  
        else:# at every round, only particle creation is recalculated
            #Parallel(n_jobs=n_proc, backend = 'loky')(delayed(RUN_FREEZE)(_) for _ in range(rounds))
            process_map(RUN_FREEZE,range(rounds),max_workers = n_proc) 
                
if __name__ == "__main__":
    sys.exit(main())        
