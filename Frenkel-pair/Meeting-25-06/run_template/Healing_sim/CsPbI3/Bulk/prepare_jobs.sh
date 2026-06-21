hm=$(pwd)
nmirror_base=500
slurm=slurm.sh
tsp_sh=tsp.sh
base="base.py"
key="FP"

rm -rf ${key}_*
nFPs=(1 10 50 100)
for FP in ${nFPs[@]}; do
	echo $FP
	f=${key}_$FP
	
	nmirror=$(echo "$nmirror_base/$FP" | bc)

	for i in $(seq 1 $nmirror);do
	path=$f/$i
	mkdir -p $path
	cp $base $path
	cp $tsp_sh $path
	cp $slurm $path
	
	cd $path
	sed -i "s|XXX|$FP|g" $tsp_sh
	sed -i "s|XXX|$FP|g" $slurm
	sbatch $slurm
	#tsp ./$tsp_sh

	cd ../../
	done
done
