import os
import shutil

def main():
	working_dir = os.getcwd()+'/'
	pip_dir = [x for x in os.popen('pip show kmc').read().split('\n') if 'Location:' in x][0].split()[1]
	pip_dir = pip_dir+'/kmc/'
	#print(working_dir)
	#print(pip_dir)
	with open(pip_dir+'pathfile.txt','w') as f:
		f.write(working_dir)
	#s = "voila Dashboard_KMC.ipynb --Voila.tornado_settings=\"{'websocket_max_message_size': 209715200}\" "
	s = "voila " +pip_dir+"Dashboard_KMC.ipynb --Voila.tornado_settings=\"{'websocket_max_message_size': 209715200}\" "
	os.system(s)
