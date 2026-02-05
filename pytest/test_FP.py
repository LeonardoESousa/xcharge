
import os
import time
import shutil

import numpy as np
import pandas as pd
import pytest

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

txts = [ x for x in os.listdir(BASE_DIR) if ".txt" in x ]
for txt in txts:
    os.remove(txt)


INPUT = f"{BASE_DIR}/FP_generation_input.py"
#CONFIGS = [10**x for x in [-5,-6,-7,-8]]
CONFIGS = [10**x for x in [-5,-6,-7]]
THRESHOLD = 1E-4

def read_dat(fn):
  data = pd.read_csv(fn,comment='#')
  data = data[data.Time != 'END']
  data.loc[(data.Time == '0'),'Time'] = 1
  return data

def average_genpop(filename):
  data = read_dat(filename)
  time = sorted(data[data['CausaMortis'] == 'generated']['Time'].to_numpy(dtype=float))
  pop = np.cumsum(np.ones([len(time)]))/10000
  x = np.mean([pop[i]/time[i] for i in range(len(pop))])
  return x
  
for conf in CONFIGS:
  os.system(f'kmc {INPUT} {conf}')


estimated_K = [ average_genpop(f'{BASE_DIR}/Simulation_generation_{conf}.txt') for conf in CONFIGS ]


rms = np.array(CONFIGS)-np.array(estimated_K)
rms = np.sqrt( sum([ x**2 for x in rms])/len(rms))

print(f'rms: {rms:.2e}')
assert rms < THRESHOLD
