import numpy as np

# Dealing with user-default options
identifier      = 'teste'
animation_mode  = False
save_animation  = False
animation_exten = "gif"
time_limit      = np.inf
pause           = False
marker_type     = 1
rotate          = False
frozen_lattice  = False
bimolec         = False       
periodic        = False
n_proc          = 1
rounds          = 1
cutoff          = np.inf
time_units      = 'ps'
generation      = {'ghost':None}
##### ANIMATION SETTINGS  ##################################
material_label = {'mat 0':0,'mat 1':1, 'mat 2': 2}
colors_dic     = {0:'black', 1:'blue', 2:'red', 3:'green', 4:'yellow'}
sizes_dic      = {0:50, 1:50, 2:50, 3:50, 4:50}
square_ratio   = True
material_leg   = True
plot_type      = "sphere" # options: sphere or scatter
scatter_alpha  = 0.3
clean_vis      = True
###############################################################


