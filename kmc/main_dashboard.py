import os
import shutil

working_dir = os.getcwd()+'/'
pip_dir = [x for x in os.popen('pip show kmc').read().split('\n') if 'Location:' in x][0].split()[1]
def main():
	print(working_dir)
	shutil.copy(pip_dir+'/kmc/Dashboard_KMC.py', working_dir)
	os.rename('Dashboard_KMC.py', 'Dashboard_KMC.ipynb')
	s = "voila Dashboard_KMC.ipynb --Voila.tornado_settings=\"{'websocket_max_message_size': 209715200}\" "
	#s = "voila Dashboard_KMC.py --Voila.tornado_settings=\"{'websocket_max_message_size': 209715200}\" "
	os.system(s)
