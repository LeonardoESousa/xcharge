import numpy as np
import random
from kmc.rates import *
from kmc.particles import *
from kmc.system import System
from kmc.bimolecular import *
from kmc.main_dashboard import main as main_dash
import sys
import warnings
import os
import copy
import matplotlib.pyplot as plt
from matplotlib import animation
import shutil
from mpl_toolkits.mplot3d import Axes3D
import importlib
warnings.filterwarnings("ignore")   
from tqdm.contrib.concurrent import thread_map, process_map
import subprocess
import inspect

#from joblib import Parallel, delayed

if sys.argv[1] == 'dash':
    main_dash()


#importing param module
working_dir = os.getcwd()+'/'
spec  = importlib.util.spec_from_file_location(sys.argv[1].split('.')[0], working_dir+sys.argv[1])
param = importlib.util.module_from_spec(spec)
spec.loader.exec_module(param)




argumentos = []  
for name, value in vars(param).items():
    if hasattr(value, 'assign_to_system') and not inspect.isclass(value): #getting only instancied classes         
        argumentos.append(value)
        #print(name,hasattr(value, 'make'))

#getting all essential info from user's input
n_proc              = param.n_proc
rounds              = param.rounds
monomolecular       = param.monomolecular
processes           = param.processes

# Dealing with user-default options
try:
    identifier     = param.identifier   
except:
    identifier = spec.name
try:
    animation_mode = param.animation_mode
except:
    animation_mode = False
try:
    save_animation = param.save_animation
except:
    save_animation = False
try:
    animation_exten = param.animation_exten
except:
    animation_exten = "gif"
try:
    time_limit = param.time_limit
except:
    time_limit = np.inf
try:
    pause = param.pause
except:
    pause = False
try:
    marker_type = param.marker_type
except:
    marker_type = 1
try:
    rotate = param.rotate
except:
    rotate = False
try:
    frozen_lattice  = param.frozen
except:
    frozen_lattice = False
try:
    bimolec  = param.bimolec
except:
    bimolec  = False       
#####


def passar(*args):
    pass

def anni(system,array,local):
    anni_general(system,array,local)
   

if bimolec:
    bi_func = anni
else:
    bi_func = passar


def regular_distance(system,local):
    dx = system.X - system.X[local]   
    dy = system.Y - system.Y[local]  
    dz = system.Z - system.Z[local]
    return dx, dy, dz

def periodic_distance(system,local):
    dx = system.X - system.X[local]   
    dy = system.Y - system.Y[local]  
    dz = system.Z - system.Z[local]
    maskx = 2*abs(dx) > system.Lx
    masky = 2*abs(dy) > system.Ly
    maskz = 2*abs(dz) > system.Lz
    dx[maskx] = -1*(system.Lx - abs(dx))[maskx]
    dy[masky] = -1*(system.Ly - abs(dy))[masky]
    dz[maskz] = -1*(system.Lz - abs(dz))[maskz]
    return dx, dy, dz

try:
    if param.periodic:
        distance = periodic_distance
    else:
        distance =  regular_distance    
except:
    distance =  regular_distance 

            
def make_system():
    #Create instance of system
    system = System()
    #Sets system properties  
    for argumento in argumentos:
        argumento.assign_to_system(system)
    system.set_basic_info(monomolecular,processes,identifier,animation_mode,time_limit,pause,bimolec,distance) 
 
    return system 
syst = make_system() #global system-type object to be used in frozen simulations, must be kept here!
    
# runs the annihilations defined in anni_funcs_array                 
def anni_general(system,anni_dict,local):
    Ss = system.particles.copy()   
    locs = np.array([s.position for s in Ss])
    superpostos = np.where(locs == local)[0]
    if len(superpostos) > 1:
        try:
            anni_dict[tuple(sorted(Ss[i].species for i in superpostos[:2]))](Ss,system,superpostos[:2])
        except:#in case of the set (particle1,partticle2) is not defined
            pass


def decision(s,system):
    kind = s.species      
    local= s.position        
    dx, dy, dz = distance(system,local)
    r = np.sqrt(dx*dx + dy*dy + dz*dz)
    hop  = system.processes[kind] 
    mono = system.monomolecular[kind]     
    jump_rate = [transfer.rate(r=r,dx=dx,dy=dy,dz=dz,system=system,particle=s) for transfer in hop]

    try:
        locais    = np.array([np.where(random.uniform(0,1) <= np.cumsum(x/np.sum(x)))[0][0] for x in jump_rate]).astype(int)
        jump_rate = np.array([jump_rate[i][locais[i]] for i in np.arange(len(locais))])
        #jump_rate = np.array([max(jump_rate[i]) for i in np.arange(len(locais))])
    except:
        locais    = np.array([local])
        jump_rate = np.array([0])

    mono_rate = np.array([m.rate(material=system.mats[local]) for m in mono])
    jump_rate = np.append(jump_rate,mono_rate)
    locais2   = np.zeros(len(mono_rate)) + local
    locais    = np.append(locais,locais2.astype(int))
    labels = hop+mono 

    jump = np.where(random.uniform(0,1) <= np.cumsum(jump_rate/np.sum(jump_rate)))[0][0]
    s.process = labels[jump]
    s.destination = locais[jump]
    return np.sum(jump_rate)

########ITERATION FUNCTIONS#######################################################
def step_ani(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        X,Y,Z = system.X, system.Y, system.Z     
        R = np.array([decision(s,system) for s in Ss])
        system.time += (1/max(R))*np.log(1/random.uniform(0,1))
        jumps = np.where(random.uniform(0,1) <= R/max(R))[0]
        for jump in jumps:
            if Ss[jump] in system.particles:
                Ss[jump].process.action(Ss[jump],system,Ss[jump].destination)   
                bi_func(system,bimolec_funcs_array,Ss[jump].destination)
        return Ss       
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1,'alive')
  
def step_nonani(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        X,Y,Z = system.X, system.Y, system.Z     
        R = np.array([decision(s,system) for s in Ss])
        system.time += (1/max(R))*np.log(1/random.uniform(0,1))
        jumps = np.where(random.uniform(0,1) <= R/max(R))[0]
        for jump in jumps:
            if Ss[jump] in system.particles:
                Ss[jump].process.action(Ss[jump],system,Ss[jump].destination)   
                bi_func(system,bimolec_funcs_array,Ss[jump].destination)       
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1,'alive')
##########################################################################################

if animation_mode:
    step = step_ani
else:
    step = step_nonani


#Prints Spectra
def spectra(system):
    if os.path.isfile("Simulation_"+system.identifier+".txt") == False:
        with open("Simulation_"+system.identifier+".txt", "w") as f:
            texto = "Time,DeltaX,DeltaY,DeltaZ,Type,Energy,Location,FinalX,FinalY,FinalZ,CausaMortis,Status"
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
        if (class_name in ["Create_Particles","Create_Particles_PROB"]):
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
