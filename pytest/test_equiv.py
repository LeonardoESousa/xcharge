#(note to myself) to run: pytest
#                 to run with print statements: pytest -s
import os
import subprocess
import pandas as pd
import numpy as np
import io
import pytest
import sys

##################################
txt_files = [i for i in os.listdir('.') if (os.path.isfile(i)) and ('.txt' in i.lower())]
for txt in txt_files:
	os.remove(txt) 
	
def get_ds(simul):
	data = pd.read_csv(simul)
	data = data[data.Time != 'END']
	dx = data['DeltaX'].to_numpy()/10
	dy = data['DeltaY'].to_numpy()/10
	dz = data['DeltaZ'].to_numpy()/10
	ds = dx*dx+dy*dy+dz*dz
	return ds

def run_job(config_test,log_name):
	lab = []
	for test in config_test:
		label = '_'.join([str(x) for x in test])
		param = ' '.join([str(x) for x in test])
		name  = 'Simulation_'+log_name.split('_')[-1].split('.')[0]+'_'+label+'.txt'
		lab.append(name)
		command = "kmc "+log_name+' '+label+' '+param
		subprocess.run(command, shell=True)
	return [(lab)]

###### EQUIVALENCE TEST ######	
log_name1    = "input_test1.py"
config_test1 = [[1000,10],[10000,1]] # 1000 rounds 10 ex, 10000 rounds 1 ex...
simuls = run_job(config_test1,log_name1)
def cond(array,cutoff):
	for x in array:
		if x > cutoff:
			return False
	return True
	
@pytest.mark.parametrize("simuls",simuls)
def test_fluor(simuls):
	avgs = []
	for simul in simuls:
		data = pd.read_csv(simul)
		data = data[data.Time != 'END']
		fluor = data['Time'].to_numpy(float)/1000
		avg = np.mean(fluor)
		avgs.append(avg)

	ref = 3
	print()
	print('test fluor')
	print(avgs)
	avgs = [ abs((av-ref)/ref) for av in avgs]
	print(avgs)
	assert cond(avgs,0.3) == True
	

@pytest.mark.parametrize("simuls",simuls)
def test_ld(simuls):
	avgs = []
	for simul in simuls:
		ds = get_ds(simul)
		ld = np.sqrt(np.mean(ds))
		avgs.append(ld)
	std = np.std(avgs)
	print()
	print('test ld')	
	print('ld: ',avgs)
	print('std:',std)
	assert std <= 0.2
	
########### RF tests #########
log_name2    = "input_test2.py"
config_test2 = [[1],[10],[25]] #1 ex, 10 ex, 100 ex
simuls2 = run_job(config_test2,log_name2)
@pytest.mark.parametrize("simuls",simuls2)
def test_Rfs(simuls):
	avgs = []
	for simul in simuls:
		ds = get_ds(simul)
		ld = np.sqrt(np.mean(ds))
		avgs.append(ld)
	std = np.std(avgs)
	print()
	print('test Rfs')
	print('ld: ',avgs)
	print('std:',std)	
	assert std <= 0.1

####### Diffusion test #######


def spectrum(dx,gran,mindx,maxdx):
	num = int((maxdx-mindx)/gran)
	if num == 0:
		bins = 1
	else:
		bins = np.linspace(mindx,maxdx,num)    
	hist, bins = np.histogram(dx,bins=bins,density=True)
	bins = bins[:-1] + (bins[1:] - bins[:-1])/2    
	return hist,bins	
def get_min_max(simuls):
	ds_arr = []
	for simul in simuls:
		ds = get_ds(simul)
		ds_arr.append(ds)
	min_ds = min([min(d) for d in ds_arr])
	max_ds = max([max(d) for d in ds_arr])
	return min_ds, max_ds
	
@pytest.mark.parametrize("simuls",simuls)
def test_diff_rounds(simuls):
	hist_arr = []
	gran = 1.5
	min_ds, max_ds = get_min_max(simuls)
	for simul in simuls:
		ds = get_ds(simul)
		hists, binss = spectrum(ds,gran,min_ds,max_ds)
		hist_arr.append(hists)
	
	ref_hist = hist_arr[0]
	N        = len(ref_hist)
	rms_arr  = []
	for i in range(1,len(hist_arr)):
		rms = np.array([(1/N)*(ref_hist[k]-hist_arr[i][k])**2 for k in range(len(ref_hist))])
		rms = np.sqrt(rms)
		rms_arr.append(rms)
	print()
	print('test diffusion rounds')
	#print('rms: ', rms)
	assert max(rms) <= 1E-3

@pytest.mark.parametrize("simuls",simuls2)
def test_diff_rf(simuls):
	hist_arr = []
	gran = 1.5
	min_ds, max_ds = get_min_max(simuls)
	for simul in simuls:
		ds = get_ds(simul)
		hists, binss = spectrum(ds,gran,min_ds,max_ds)
		hist_arr.append(hists)
	
	ref_hist = hist_arr[0]
	N        = len(ref_hist)
	rms_arr  = []
	for i in range(1,len(hist_arr)):
		rms = np.array([(1/N)*(ref_hist[k]-hist_arr[i][k])**2 for k in range(len(ref_hist))])
		rms = np.sqrt(rms)
		rms_arr.append(rms)
	print()
	print('test diffusion Rf')
	#print('rms: ', rms)
	assert max(rms) <= 1E-3

