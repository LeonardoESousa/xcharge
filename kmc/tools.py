# Aux. code for analysis 
import numpy as np
import ipywidgets as widgets

##############################
### FUNCTIONS FOR DASHBOARD ##
##############################
#sets the bounds for some variable between [min_val,max_val]
def get_bound(max_val,min_val):
    if max_val  < 0:
        boundary_max   = -.1 
    else:
        boundary_max   = .1
    if  min_val  < 0:
        boundary_min   = .1
    else:
        boundary_min   = -.1
    return min_val*(1+boundary_min),max_val*(1+boundary_max)

##################################################################

#filter some variable inside some list
def filter_data_list(data,keyword,list_choosen):
    return data[(data[keyword].isin(list_choosen))]

#filter some variable inside some numeric range
def filter_data_range(data,keyword,values):
    if keyword != 'Time':
        return data[((data[keyword].to_numpy(float) <= values[1]) & (data[keyword].to_numpy(float) >= values[0]))] 
    else:
        return data[((data[keyword].to_numpy(float) <= 10**(12+values[1])) & (data[keyword].to_numpy(float) >= 10**(12+values[0])))]

#Returns a filtered dataframe for a given particle
#death processes, material location, energy...
def filter_per_particle(dataframe,particle,deaths,mats,energy,times):
    data = filter_data_list(dataframe,'Type',[particle])
    data_filtred = data
    data_filtred = filter_data_range(data_filtred,'Time',times)
    data_filtred = filter_data_list(data_filtred,'CausaMortis',deaths)
    data_filtred = filter_data_range(data_filtred,'Location',mats)
    data_filtred = filter_data_range(data_filtred,'Energy',energy)
    #data_filtred = filter_data_list(data_filtred,'Status',status)   
    return data_filtred
###################################################################


#returns a set of widgets for some particle to compose
#the first panel
def get_widgets_per_particle(min_energy,max_energy,unique_mats,unique_death,min_time,max_time):
    min_bound,max_bound = get_bound(max_energy,min_energy)
    energy = widgets.FloatRangeSlider(
        value=[min_energy, max_energy],
            min=min_bound,
            max=max_bound,
            step=(max_energy-min_energy)/4,
            description='Energy (eV):',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='.1f',
        )       
    mats = widgets.IntRangeSlider(
            value=[np.amin(unique_mats), np.amax(unique_mats)],
            min=np.amin(unique_mats),
            max=np.amax(unique_mats),
            step=1,
            description='Materials:',
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='d',
        )                
    deaths = widgets.SelectMultiple(
            options=unique_death,
            value=[unique_death[0]],
            #rows=10,
            description='Process',
            disabled=False,
        )
    #status = widgets.SelectMultiple(
    #        options=unique_status,
    #        value=[unique_status[0]],
    #        #rows=10,
    #        description='Status',
    #        disabled=False,
    #    )
    times = widgets.FloatRangeSlider(
            value=[np.log10(min_time)-12, np.log10(max_time)-12],
            min=np.log10(min_time)-12,
            max=np.log10(max_time)-12,
            step=0.1,
            description='Time: $10^x$(s)',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='.1f',
        )
    return [mats,energy,deaths,times]#,status]


#gathers the personalized info to construct the widgets
#from function get_widgets_per_particle and returns it
def join_widgets(dataframe,particle):
    
    data = filter_data_list(dataframe,'Type',[particle])
    unique_mats   = data['Location'].unique()
    unique_death  = data['CausaMortis'].unique()
    unique_energy = data['Energy'].unique()
    unique_status = data['Status'].unique()
    max_energy = np.amax(unique_energy)
    min_energy = np.amin(unique_energy)
    energy, mats, deaths,status = get_widgets_per_particle(min_energy,max_energy,
                                                           unique_mats,unique_death,unique_status)
    return [particle,energy, mats, deaths,status]
    
#wraps death,status,mats and enegy widgets in a single big box  
def make_boxes(wids):
    deaths,mats,energy,times = wids[2],wids[0],wids[1], wids[3],
    left_box  = widgets.VBox([deaths])#,status])
    right_box = widgets.VBox([mats,energy,times])
    ui = widgets.HBox([left_box,right_box])  
    return ui    
    
def trpl(times,bin_num=10):
    num = int((max(times)-min(times))/bin_num)
    hist, bins = np.histogram(times, bins=np.linspace(min(times),max(times),num),density=True)    
    bins = bins[:-1] +(bins[1:] - bins[:-1])/2
    return  hist, bins

    
def spectrum(dx,gran):
    num = int((max(dx)-min(dx))/gran)
    if num == 0:
        bins = 1
    else:
        bins = np.linspace(min(dx),max(dx),num)    
    hist, bins = np.histogram(dx,bins=bins,density=True)
    bins = bins[:-1] + (bins[1:] - bins[:-1])/2    
    return hist,bins
    
def drift(data):
    t  = data['Time'].to_numpy(dtype=float)
    dx = data['DeltaX'].to_numpy()
    dy = data['DeltaY'].to_numpy()
    dz = data['DeltaZ'].to_numpy()
    mux = np.mean(dx/t)
    muy = np.mean(dy/t)
    muz = np.mean(dz/t)
    return np.array([mux,muy,muz])

def get_peak(hist,bins):
    ind   = np.where(hist == max(hist))[0][0]
    peak  = bins[ind]
    return peak



################  
# HISTOGRAM FUNCS

#generate N_boxes for a given particle (Histogram functionallity)    
def get_boxes_for_a_particle(N_boxes,death_list_particle):
    death_list_particle = [x for x in death_list_particle if x != 'Absent']
    n_deaths = len(death_list_particle)
    real_boxes_indx = []
    boxes = []
    for i in range(N_boxes):
        if i+1 <= n_deaths:
            death_name = death_list_particle[i]
            boxes.append(widgets.Checkbox(value=False,description=str(death_name),disabled=False))
        else:
            boxes.append(widgets.Checkbox(value=False,description='---',disabled=True))
    return boxes,death_list_particle    
    
    
    
    
    
    
   

