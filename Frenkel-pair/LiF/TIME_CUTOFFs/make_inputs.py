#######################
import os
import numpy as np
import shutil
from itertools import product
#####################
hm     = os.getcwd()
base   = "base.py"
cif    = "LiF.cif"
vals   = [0.1, 0.5, 0.7, 1.0]#[0.1, 0.5, 0.7, 1.0]
orders = [1E-4,1E-3]#[1E-4,1E-3,1E-2]
times  = [ f"{x[0]*x[1]:.2e}" for x in  product(vals, orders)]
times  = list(set(times))

keyword="TIMES_"
try:
    dirs   = [x for x in os.listdir() if keyword in x ]
    [ shutil.rmtree(x) for x in dirs]
except Exception as e:
    print(e)
    pass
for time in times:
    fname= f"{keyword}{time}"
    os.mkdir(fname)
    shutil.copy(f"{hm}/{base}",fname)
    shutil.copy(f"{hm}/{cif}",fname)
    os.chdir(fname)
    command = f"ts kmc base.py {time}"
    os.system(command)
    os.chdir('..')

