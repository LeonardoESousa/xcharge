import numpy as np
import random
from kmc.system import System
import kmc.bimolecular
import sys
import warnings
warnings.filterwarnings("ignore") 
import os
import copy
import matplotlib.pyplot as plt
from matplotlib import animation
import importlib  
from tqdm.contrib.concurrent import thread_map, process_map
import inspect
import kmc.variables

#from joblib import Parallel, delayed

#importing param module
spec  = importlib.util.spec_from_file_location(sys.argv[1].split('.')[0], os.path.join(os.getcwd(),sys.argv[1]))
param = importlib.util.module_from_spec(spec)
spec.loader.exec_module(param)




argumentos = []  
for name, value in vars(param).items():
    if hasattr(value, 'assign_to_system') and not inspect.isclass(value): #getting only instancied classes         
        argumentos.append(value)
        #print(name,hasattr(value, 'make'))

#getting all essential info from user's input
monomolecular       = param.monomolecular
processes           = param.processes

def set_variables(name):
    try:
        return getattr(param, name)
    except:
        return getattr(kmc.variables, name)

identifier     = set_variables('identifier')     
animation_mode = set_variables('animation_mode') 
save_animation = set_variables('save_animation') 
animation_exten= set_variables('animation_exten')
time_limit     = set_variables('time_limit')     
pause          = set_variables('pause')          
marker_type    = set_variables('marker_type')    
rotate         = set_variables('rotate')         
frozen_lattice = set_variables('frozen_lattice') 
bimolec        = set_variables('bimolec')  
periodic       = set_variables('periodic')          
n_proc         = set_variables('n_proc')
rounds         = set_variables('rounds')
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
    dx[maskx] = -1*np.sign(dx[maskx])*(system.Lx - abs(dx))[maskx]
    dy[masky] = -1*np.sign(dy[masky])*(system.Ly - abs(dy))[masky]
    dz[maskz] = -1*np.sign(dz[maskz])*(system.Lz - abs(dz))[maskz]
    return dx, dy, dz

if periodic:
    distance = periodic_distance
else:
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
    mono_rate = [m.rate(material=system.mats[local]) for m in mono]
    jump_rate.append(mono_rate)
    sizes     = np.array([len(i) for i in jump_rate])
    jump_rate = np.concatenate(jump_rate)
    labels    = hop+mono 
    soma      = np.sum(jump_rate)
    jump      = np.argmax(random.uniform(0,1) <= np.cumsum(jump_rate/soma))
    destino   = np.argmax(np.cumsum(sizes)-1 >= jump)
    s.process = labels[destino]
    if destino < len(hop):
        s.destination = jump - np.sum(sizes[:destino])
    else:
        s.destination = local
    return soma

########ITERATION FUNCTIONS#######################################################
''' 
def step_ani(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        R = np.array([decision(s,system) for s in Ss])
        system.time += np.mean((1/R)*np.log(1/random.uniform(0,1)))
        for s in Ss:
            if s in system.particles:
                s.process.action(s,system,s.destination)   
                bi_func(system,kmc.bimolecular.bimolec_funcs_array,s.destination)
        return Ss       
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1,'alive')
 
def step_nonani(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        R = np.array([decision(s,system) for s in Ss])
        system.time += np.mean((1/R)*np.log(1/random.uniform(0,1)))
        #ideia annirad: olhar cada par de particulas e ver se o raio corresponde ao que foi colocado no input, se não, arrumar
        for s in Ss:
            if s in system.particles:
                s.process.action(s,system,s.destination)   
                bi_func(system,kmc.bimolecular.bimolec_funcs_array,s.destination)
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1,'alive')
'''
'''
#ideia        
def step_nonani(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        #ideia annirad: olhar cada par de particulas e ver se o raio corresponde ao que foi colocado no input, se não, arrumar
        R = []
        for s in Ss:
            if s in system.particles:
                Rs = decision(s,system)
                checar se tem particula onde vou pular
                se sim, voltar pro Rs. Se n continua
                s.process.action(s,system,s.destination)
                bi_func(system,kmc.bimolecular.bimolec_funcs_array,s.destination)
    system.time += np.mean((1/R)*np.log(1/random.uniform(0,1)))   
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1,'alive')
'''      
def step_nonani(system): 
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        R = []
        for s in Ss:
            if s in system.particles:
                Rs = decision(s,system)
                s.process.action(s,system,s.destination)
                bi_func(system,kmc.bimolecular.bimolec_funcs_array,s.destination)
                R.append(Rs)
        R = np.array(R)
        system.time += np.mean((1/R)*np.log(1/random.uniform(0,1)))   
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1,'alive')
 
def step_ani(system):
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        R = []
        for s in Ss:
            if s in system.particles:
                Rs = decision(s,system)
                s.process.action(s,system,s.destination)
                bi_func(system,kmc.bimolecular.bimolec_funcs_array,s.destination)
                R.append(Rs)
        R = np.array(R)
        system.time += np.mean((1/R)*np.log(1/random.uniform(0,1)))
        return Ss
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
    filename = f"Simulation_{system.identifier}.txt"
    if os.path.isfile(filename) == False:
        with open(filename, "w") as f:
            texto = "Time,DeltaX,DeltaY,DeltaZ,Type,Energy,Location,FinalX,FinalY,FinalZ,CausaMortis,Status"
            f.write(texto+"\n") 
    texto = ''
    for s in system.dead:
        texto += s.write()
    with open(filename, "a") as f:   
        f.write(texto+'END\n')
        
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
                                    interval=150, blit=False,repeat=False,cache_frame_data=True)#,save_count=1000)  
                             
                                       
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
