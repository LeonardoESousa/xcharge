import numpy as np
import pandas as pd

data = pd.read_csv('Simulation_forster_singlet.txt', delim_whitespace=True)
data = data[data.Time != 'END'] #removes the end of round lines
data = data[data.Location == 0]
tempos = np.array(data['Time'],dtype='float')
print(np.mean(tempos))

