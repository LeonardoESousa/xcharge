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
from matplotlib.lines import Line2D
import importlib  
import inspect
import kmc.variables
import kmc.utils
try:
    from importlib import metadata
except ImportError: # for Python<3.8
    import importlib_metadata as metadata
import multiprocessing 
import tqdm
# Setting up the interface
print('####################################################################')
print('Xcharge: A Kinetic Monte Carlo Model for Exciton and Charge Dynamics')
#print('Repo link: '+metadata.metadata('kmc')['Home-page'])
#print('Authors  : '+ metadata.metadata('kmc')['Author'])
print('Version  : '+  metadata.metadata('kmc')['VERSION'])
print()
output_header='# Version  : '+  metadata.metadata('kmc')['VERSION']+'\n'
### end interface


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

generation     = set_variables('generation')
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
cutoff         = set_variables('cutoff')
FPcutoff       = set_variables('FPcutoff')
colors_dic     = set_variables('colors_dic')
sizes_dic      = set_variables('sizes_dic')
square_ratio   = set_variables('square_ratio')
scatter_alpha  = set_variables('scatter_alpha')
clean_vis      = set_variables('clean_vis')
material_label = set_variables('material_label')
material_leg   = set_variables('material_leg')
sizes_dic      = set_variables('sizes_dic')
plot_type      = set_variables('plot_type')
time_units     = set_variables('time_units')
particle_condition = set_variables('particle_condition')
print_site_position = set_variables('print_site_position')
#####


def passar(*args):
    pass

def anni(system,array,local):
    anni_general(system,array,local)
   

if bimolec:
    bi_func = anni
else:
    bi_func = passar


def regular_distance(system,local,destination=None):
    if destination is not None:
        dx = system.X[destination] - system.X[local]   
        dy = system.Y[destination] - system.Y[local]  
        dz = system.Z[destination] - system.Z[local]
    else:
        #dx,dy,dz = kmc.utils.distance(system.X,system.Y,system.Z,len(system.X),local)
        dx = system.X - system.X[local]   
        dy = system.Y - system.Y[local]  
        dz = system.Z - system.Z[local]
    return dx, dy, dz

def periodic_distance(system,local,destination=None):
    if destination is None:
        dx = system.X - system.X[local]   
        dy = system.Y - system.Y[local]  
        dz = system.Z - system.Z[local]
    else:
        dx = system.X[destination] - system.X[local]   
        dy = system.Y[destination] - system.Y[local]  
        dz = system.Z[destination] - system.Z[local]
    if system.Lx > 0:
        dx -= system.Lx*np.round(dx/system.Lx)   
    if system.Ly > 0:
        dy -= system.Ly*np.round(dy/system.Ly)  
    if system.Lz > 0:
        dz -= system.Lz*np.round(dz/system.Lz)
    return dx,dy,dz

if periodic:
    distance = periodic_distance
else:
    distance =  regular_distance    

            
def make_system():
    #Create instance of system
    system = System()
    #Sets system properties  
    system.set_basic_info(monomolecular,processes,identifier,animation_mode,time_limit,pause,bimolec,distance,generation,[cutoff,FPcutoff]) 
    for argumento in argumentos:
        argumento.assign_to_system(system)
 
    return system 
    
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
    r = kmc.utils.distances(dx,dy,dz,len(dx))
    cut = np.where(r < cutoff)[0]
    r = r[cut]
    dx = dx[cut]
    dy = dy[cut]
    dz = dz[cut]
    mats = system.mats[cut]
    hop  = system.processes[kind] 
    mono = system.monomolecular[kind]     
    jump_rate = [transfer.rate(r=r,dx=dx,dy=dy,dz=dz,system=system,particle=s,mats=mats,matlocal=system.mats[local],cut=cut) for transfer in hop] 
    mono_rate = [[m.rate(material=system.mats[local])] for m in mono]
    jump_rate.extend(mono_rate)   
    sizes     = np.array([len(i) for i in jump_rate])
    jump_rate = np.concatenate(jump_rate)
    labels    = hop+mono
    soma, jump = kmc.utils.jump(np.array(jump_rate,dtype=np.float64),len(jump_rate),random.uniform(0,1)) #dtype mismatch if the number, by chance, is integer. Force conversion
    destino   = np.argmax(np.cumsum(sizes) -1 >= jump)
    s.process = labels[destino]
    if destino < len(hop):
        s.destination = cut[int(jump - np.sum(sizes[:destino]))]
        if s.process.__class__.__name__ == "Formation": #different cut, different handling
            #fix here needed, failed attempt below
            #FPCUT=np.where(r < system.FPcutoff)[0]
            #s.destination = FPCUT[int(jump - np.sum(sizes[:destino]))]
            pass
    else:
        s.destination = local
    return soma

########ITERATION FUNCTIONS#######################################################
def chose_generation(system):
    try:
      genrates = system.generation
      K_gen = [ gen.rate() for gen in genrates ]
      prob_gens   = np.cumsum(K_gen)/np.sum(K_gen)
      u = random.uniform(0,1)
      choose_gen = np.argmax(prob_gens >= u)
      return genrates[choose_gen], K_gen[choose_gen]
    except Exception as e:
      print("Generation errror:",e)
      return 0,0
def step_nonani(system):
    while (((not particle_condition) or (system.count_particles() > 0)) and (system.time < system.time_limit)): # if particle_condition = True  this equals to system.count_particles() > 0 and system.time < system.time_limit
        #print(system.IT,f'{system.time:.2e}',len(system.particles))
        system.IT += 1
        event_g, K_g = chose_generation(system)
        
        Ss = system.particles.copy()
        Rs = np.array([decision(s, system) for s in Ss])
        sum_Rs = np.sum(Rs)
        
        R_total = sum_Rs + K_g
        Prob    = [K_g,sum_Rs]
        Prob    = np.cumsum(Prob)/R_total
        
        dt = (1 / R_total) * np.log(1 / random.uniform(0, 1))
        system.time += dt

        #print(Prob)
        #print(f'IT {system.IT} before:',[[s.species,s.status,s.position,s.ghost_site] for s in system.particles if s.species == "frenkelpair"])
        u = random.uniform(0, 1)
        if u < Prob[0]:
            # --- GENERATION EVENT ---
            event_g.action(system, npart=1)
            #print(f'at time {system.time:.10e}, event: generation, dt: {dt:.2e}')
        else:
            # --- PARTICLE-BASED EVENT ---
            u_part = random.uniform(0, 1)#u - K_g 
            prob_part = np.cumsum(Rs)/sum_Rs
            choose_part = np.argmax(prob_part >= u_part)
            
            s = Ss[choose_part]
            
            # Execute the action
            s.process.action(s, system, s.destination)
            bi_func(system, kmc.bimolecular.bimolec_funcs_array, s.destination)
            s.stamp_time(system)
            #print(prob_part,u,prob_part[choose_part])
            #print(f'at time {system.time:.10e}, event: {s.process}, dt: {dt:.2e}')
            #print([system.mats[s.position] for s in system.particles])
        #print(f'IT {system.IT} after:',[[s.species,s.status,s.position,s.ghost_site] for s in system.particles if s.species == "frenkelpair"])             
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1,'alive')
        s.stamp_time(system)

def step_ani(system):
    while (((not particle_condition) or (system.count_particles() > 0)) and (system.time < system.time_limit)): # if particle_condition = True  this equals to system.count_particles() > 0 and system.time < system.time_limit
        #print(system.IT,f'{system.time:.2e}',len(system.particles))
        system.IT += 1
        event_g, K_g = chose_generation(system)
        
        Ss = system.particles.copy()
        Rs = np.array([decision(s, system) for s in Ss])
        sum_Rs = np.sum(Rs)
        
        R_total = sum_Rs + K_g
        Prob    = [K_g,sum_Rs]
        Prob    = np.cumsum(Prob)/R_total
        
        dt = (1 / R_total) * np.log(1 / random.uniform(0, 1))
        system.time += dt

        #print(Prob)
        #print(f'IT {system.IT} before:',[[s.species,s.status,s.position,s.ghost_site] for s in system.particles if s.species == "frenkelpair"])
        u = random.uniform(0, 1)
        if u < Prob[0]:
            # --- GENERATION EVENT ---
            event_g.action(system, npart=1)
            #print(f'at time {system.time:.10e}, event: generation, dt: {dt:.2e}')
        else:
            # --- PARTICLE-BASED EVENT ---
            u_part = random.uniform(0, 1)#u - K_g 
            prob_part = np.cumsum(Rs)/sum_Rs
            choose_part = np.argmax(prob_part >= u_part)
            
            s = Ss[choose_part]
            
            # Execute the action
            s.process.action(s, system, s.destination)
            bi_func(system, kmc.bimolecular.bimolec_funcs_array, s.destination)
            s.stamp_time(system)
            #print(prob_part,u,prob_part[choose_part])
            #print(f'at time {system.time:.10e}, event: {s.process}, dt: {dt:.2e}')
            #print([system.mats[s.position] for s in system.particles])
        #print(f'IT {system.IT} after:',[[s.species,s.status,s.position,s.ghost_site] for s in system.particles if s.species == "frenkelpair"])
        return system.particles.copy() #Ss  <--- returns wrong number if generation is active                
    Ss = system.particles.copy()
    for s in Ss:
        s.kill('alive',system,system.s1,'alive')
        s.stamp_time(system)
##########################################################################################

if animation_mode:
    step = step_ani
else:
    step = step_nonani

def open_log():
    filename = f"Simulation_{identifier}.txt"
    if os.path.isfile(filename) == False:
        with open(filename, "w") as f:
            f.write(output_header)
            texto = "Time,DeltaX,DeltaY,DeltaZ,Type,Energy,Location,FinalX,FinalY,FinalZ,CausaMortis,Status"
            f.write(texto+"\n") 
    return filename

#Prints Spectra
def spectra(system,f):
    texto = ''
    for s in system.dead:
        texto += s.write()
    f.write(texto+f'END\n')


def draw_sphere(ax, center, radius, color, margin_size, resolution=30):
    [[x_min,y_min,z_min],[x_max,y_max,z_max]] = margin_size


    u = np.linspace(0, 2*np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)

    x = (radius*(x_max-x_min))* np.outer(np.cos(u), np.sin(v)) + center[0]
    y = (radius*(y_max-y_min)) * np.outer(np.sin(u), np.sin(v)) + center[1]
    z = (radius*(z_max-z_min)) * np.outer(np.ones_like(u), np.cos(v)) + center[2]

    ax.plot_surface(x, y, z, color=color, shade=True)


def draw_frankelpair(rS,rG,s,ax,marker,color,size,alpha):
    xs,ys,zs = rS
    xg,yg,zg = rG
    ax.plot([xs,xg],[ys,yg],[zs,zg],color=s.color,alpha=0.5,label=s.species,linewidth=10)
def animate(num,system,ax,marker_option,rotate,colors_dic,margin_size):
    #try: #this is to freeze the frame at each iteration. Good for debug. Keeping here!
    #   animate.event_source.stop()
    #except:
    #   pass
    Ss = step(system)
    if Ss is None:
       Ss = []
    
    X, Y, Z = system.X, system.Y, system.Z        
    mats = system.mats                            
    ax.clear()
    if print_site_position:
      [ax.text(X[i],Y[i],Z[i],i) for i in range(len(X))]
    #ploting the sites according to mat index
    n_mats = [ x for x in np.unique(mats) if x != 999]
    for mat in n_mats:
        X_mat = X[mats == mat]
        Y_mat = Y[mats == mat]
        Z_mat = Z[mats == mat]
        ax.set_proj_type('persp')
        
        if plot_type == "sphere":
            indx = len(X_mat)
            for i in range(indx):
                #draw_sphere(ax, [X_mat[i],Y_mat[i],Z_mat[i]], 0.03*(sizes_dic[int(mat)]/np.amax(list(sizes_dic.values()))), colors_dic.get(int(mat)), margin_size, resolution=40)
                draw_sphere(ax, [X_mat[i],Y_mat[i],Z_mat[i]], 0.03, colors_dic.get(int(mat)), margin_size, resolution=15)
        else:
            ax.scatter(X_mat,Y_mat,Z_mat,alpha=scatter_alpha,color=colors_dic.get(int(mat)),s=sizes_dic[int(mat)],depthshade=True)
    '''#keeping this commented for a while so I can debug properly
    try:  
        for s in Ss:
            xs = X[s.position]        	
            ys = Y[s.position]
            zs = Z[s.position]    
            if s.species == 'frenkelpair':
                xg = X[s.ghost_site]        	
                yg = Y[s.ghost_site]
                zg = Z[s.ghost_site]
                ax.scatter(xg,yg,zg,alpha=scatter_alpha,color=colors_dic.get(int(s.origin_site)),s=sizes_dic[int(s.origin_site)],depthshade=True) #drawing deactivated site
                draw_frankelpair([xs,ys,zs],[xg,yg,zg],s,ax,marker=s.marker,color=s.color,size=200,alpha=1) #drawing the pair
                continue
            #print(system.mat[s.position])                
            if marker_option == 1:
                ax.scatter(xs,ys,zs,marker=s.marker,color=s.color,s=200,alpha=1,label=s.species)      
            if marker_option == 0:
                ax.scatter(xs,ys,zs,color=s.color,s=100,alpha=1,label=s.species)                 
    except Exception as e:
        print('drawing error:',e)
        pass
    '''
    ########################
    Ss =[ s for s in Ss if s.status=='alive']
    for s in Ss:
        xs = X[s.position]        	
        ys = Y[s.position]
        zs = Z[s.position]    
        if s.species == 'frenkelpair':
            xg = X[s.ghost_site]        	
            yg = Y[s.ghost_site]
            zg = Z[s.ghost_site]
            #ax.scatter(xg,yg,zg,alpha=scatter_alpha,color=colors_dic.get(int(s.origin_site)),s=sizes_dic[int(s.origin_site)],depthshade=True) #drawing deactivated site
            draw_frankelpair([xs,ys,zs],[xg,yg,zg],s,ax,marker=s.marker,color=s.color,size=200,alpha=1) #drawing the pair
            continue
            #print(system.mat[s.position])                
        if marker_option == 1:
            ax.scatter(xs,ys,zs,marker=s.marker,color=s.color,s=200,alpha=1,label=s.species)      
        if marker_option == 0:
            ax.scatter(xs,ys,zs,color=s.color,s=100,alpha=1,label=s.species)   
    ############################              
    if rotate:#rotating the animation by an angle of IT
        ax.view_init(azim = system.IT)
    #removing duplicates on legend    
    handles, labels = plt.gca().get_legend_handles_labels()
    # De-duplicate first (dict keeps last handle per label), then sort by label.
    # Sorting the raw zip breaks when labels repeat because Python tries to
    # compare the (non-orderable) artist objects to break ties.
    by_label = dict(zip(labels, handles))
    by_label = dict(sorted(by_label.items(), key=lambda kv: kv[0]))
    particle_legend = ax.legend(by_label.values(), by_label.keys())
    ax.add_artist(particle_legend)
    #ax.text2D(0.03, 0.97, "time = %.2e ps" % (system.time), transform=ax.transAxes) #time
    ax.text2D(0.03, 0.97, f"time = {system.time:.2e} {time_units}", transform=ax.transAxes) #time
    ax.text2D(0.03, 0.93, "npart  = %.0f"  % (len(system.particles)), transform=ax.transAxes) #npart
             
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')     
    
    if square_ratio:
        ax.set_box_aspect([1, 1, 1])
    if clean_vis: #cleaning the visualization
        ax.set_axis_off()
        ax.grid(False)
        ax.set_position([0, 0, 1, 1]) #animation occupies more of the screen
    if material_leg:#here, we add a specific legend for the material species, which is a fixed set of elements
        legend_elements = []
        mat_lab = material_label.keys()
        for mat in mat_lab:
          leg = Line2D([0], [0],
               marker='o',
               linestyle='None',
               label=mat,
               markerfacecolor=colors_dic[material_label[mat]],
               markeredgecolor=colors_dic[material_label[mat]],
               markersize=10)
          legend_elements.append(leg)
        mat_legend = ax.legend(handles=legend_elements, title='Materials', loc='lower right')
        ax.add_artist(mat_legend)
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
        if (class_name in ["Create_Particles","Create_Particles_PROB","Create_ParticlesFP"]):
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
    return system

def RUN_FREEZE(dynamic): #ROUND DYNAMICS WHERE, FOR EACH INSTANCE, THE LATTICE REMAINS INTACT
    syst = dynamic[1]
    system = reroll_system(copy.deepcopy(syst)) 
    step(system)
    return system


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
    
    
    #calculating the border for visualization purposes

    p_max = [ np.amax(x) for x in [ list(system.X), list(system.Y),list(system.Z)]]
    p_min = [ np.amin(x) for x in [ list(system.X), list(system.Y),list(system.Z)]]
                    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fig.canvas.mpl_connect('button_press_event', onClick) #pausing if clicking
    fig.canvas.mpl_connect('draw_event', lambda event: pause_plot(event, pause)) #pausing if pause = True at the first frame
       
    #ani = animation.FuncAnimation(fig, animate, fargs=[system,ax,marker_type,rotate,colors_dic,[p_min,p_max]],
    #                                interval=25, blit=True,repeat=False,cache_frame_data=True)#,save_count=1000)  
                             
    ani = animation.FuncAnimation(fig, animate, fargs=[system,ax,marker_type,rotate,colors_dic,[p_min,p_max]],
                                    interval=1, blit=False,repeat=False,cache_frame_data=True,save_count=1000)
                                    
    #note for future improvement: blit = true -> good for save animation, blit=false -> good for realtime interaction                                   
    animate.event_source = ani.event_source
    return ani 
    

    
def main():
        
    if animation_mode:
        ani = run_animation()
        path=identifier+"_animation."+animation_exten
                                                   
        if save_animation:                   
            
            #save .gif
            if animation_exten == 'gif':
                ani.save(path, writer='pillow', fps=10)
            
            #save .mp4
            if animation_exten == 'mp4':
                writervideo = animation.FFMpegWriter(fps=10,extra_args=["-crf", "10","-preset", "slow","-pix_fmt", "yuv420p"]) 
                ani.save(path, writer=writervideo)
        
        plt.show()
    else:      
        p = multiprocessing.Pool(n_proc) 
        filename = open_log()      
        if not frozen_lattice: # at every round, the entire lattice is recalculated
            run = RUN
            args = [(i) for i in range(rounds)]
        else:# at every round, only particle creation is recalculated
            syst = make_system()
            run = RUN_FREEZE
            args = [(i, syst) for i in range(rounds)]
        #debug
        #'''
        #for arg in args:
        #    result = run(arg)
        #    with open(filename, "a") as f:
        #        spectra(result,f)
        #'''
        with open(filename, "a") as f:
            for result in tqdm.tqdm(p.imap(run, args),total=rounds):
                spectra(result,f)        
        #'''

if __name__ == "__main__":
    sys.exit(main())        
