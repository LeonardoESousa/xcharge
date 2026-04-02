
import os
import time
import shutil

import numpy as np
import pandas as pd
import pytest

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

#txts = [ x for x in os.listdir(BASE_DIR) if ".txt" in x ]
#for txt in txts:
#    os.remove(txt)

#analytical values
W          = 2E+5 # 10.1103/PhysRevLett.57.1781 bulk diffusion rate (W) = 2*10^5 Hz (699 K)
dx         = (2.887E-8) # radius of the first neigh
coord      = 12 # coordination for this dx
dt         = (1/(coord*W))
D_analytic = (1/6)*(dx**2/dt)
############################


INPUT     = f"{BASE_DIR}/LiF_input.py"
#CONFIGS  = [10**x for x in [-5,-6,-7,-8]]
THRESHOLD = 5E-2 # 5% error on D

def read_dat(fn):
  data = pd.read_csv(fn,comment='#')
  data = data[data.Time != 'END']
  data.loc[(data.Time == '0'),'Time'] = 1
  return data

def calc_D(list_particle):
  D = np.zeros([len(list_particle),3])
  for i,(time,dx,dy,dz) in enumerate(list_particle):
      #print(i,time,dx,dy,dz)
      D[i] = [dx**2/time,dy**2/time,dz**2/time]
  D_m = (1/6)*np.sum(D,axis=1) #this is in AA^2/s
  D_m = D_m*(1E-16) #converting AA to cm
  D_m = np.mean(D_m)
  return D_m
  
def diffusion_anal(filename):
  data = read_dat(filename)
  DIFF_SLICE = ['Time','DeltaX','DeltaY','DeltaZ']
  interstitial = data[(data['CausaMortis'] == 'alive') & (data['Type'] == 'interstitial')][DIFF_SLICE].to_numpy(dtype=float)
  D_m = calc_D(interstitial)
  return D_m

os.system(f'kmc {INPUT}')

print("run finished")
print(os.listdir())

estimated_D = diffusion_anal(f'Simulation_LiF.txt')

D_diff = abs(estimated_D-D_analytic)

print(f'D_analytic = {D_analytic}\nD_calc= {estimated_D:.2e}')
print(f'D_analytic - D_calc: {D_diff:.2e}\n D diff %{D_diff/D_analytic:.2e}')
#assert D_diff < THRESHOLD

def test_D():
    assert (D_diff/D_analytic) < THRESHOLD
