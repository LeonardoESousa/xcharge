import os
import subprocess
import sys

#importing param
working_dir = os.getcwd()
sys.path.insert(1, working_dir)
import PARAM

animation_mode = PARAM.animation_mode
rounds = PARAM.rounds
n_proc = PARAM.n_proc
n_loops = int(rounds/n_proc)
n_loops_leftover= rounds % n_proc

bib_path = "XXX"
main_script = bib_path +"newkmc.py"

def run_parallel(n_jobs):
	if n_jobs > 0:
		for i in range(n_jobs):
			p = subprocess.Popen(['python3', main_script, str(working_dir)])
		p.communicate()
	else:
		pass
	
	#while p.poll() is None:
		#print('Still sleeping')
	#	time.sleep(10)
if animation_mode:	 	
	subprocess.call(['python3', main_script, str(working_dir)])	
else:
	for i in range(n_loops):	    
		run_parallel(n_proc)
		print(i,"entrou!")
	run_parallel(n_loops_leftover)
