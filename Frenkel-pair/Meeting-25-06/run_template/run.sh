hm=$(pwd)

sh1=prepare_jobs.sh
tspsh=tsp.sh
slurm=slurm.sh

cd Healing_sim
for specie in $(ls -vd */); do
	cd $specie
	for lattice in $(ls -vd */);do
	cd $lattice	
	pwd
	
	cp ${hm}/${sh1} .
	cp ${hm}/${tspsh} .
	cp ${hm}/${slurm} .

	./${sh1}

	cd ../
	done
	cd ../
done
cd ../

