#!/bin/bash
# Script that installs KMC software

hm=$(pwd)

read -s -p "Enter Password for sudo: " sudoPW

echo "--------------------------------------------------------------------------------------------------------------"
echo "Install script for KMC software LCC-UnB (2021)"
echo "Usage: Linux (Local/Global)"
echo "In addition to this script, make sure that your python has the following libs:"
echo "numpy, matplotlib, collections"
echo "It also requires task-spooler installed"
echo "--------------------------------------------------------------------------------------------------------------"
echo
echo "You can install this software by two procedures: 1) Local install 2) Installing a global package. Which one do you wish to try?"
read ans

source_dir="$hm/source" 
param_file="PARAM.py"
paralell_file="paralellizer.py"


param_bib=$(grep "bib_path =" $param_file)
paralell_bib=$(grep "bib_path =" "$source_dir/$paralell_file")





if [ $ans = "1" ]; then
	echo $ans
	sed  -i "s|$param_bib|bib_path = \"$source_dir/\"|g" $param_file
	sed  -i "s|$paralell_bib|bib_path = \"$source_dir/\"|g" "$source_dir/$paralell_file"

	echo
	echo "Add an alias in .bashrc? (y/n)"
	read bash_ans
	if [ $bash_ans = "y" ]; then
	
		echo	
		echo "Now we are going to edit your .bashrc file. Please insert your username (ex.: padawan, john..)"
		read username
		bash_path="/home/$username"

                #package_path="/home/$username/KMC_source"
                package_path="/home/tiago/Desktop/pack"
                mkdir -p $package_path
                mv -f $source_dir/$paralell_file $package_path

		if grep -q "alias kmc" "/home/$username/.bashrc"
		then
		
                        echo "I found the following entry in .bashrc with kmc keyword:"
                       	echo
		       	previous_kmc=$(grep "alias kmc" "$bash_path/.bashrc")
                        echo $previous_kmc
                        echo
			echo "Overwriting..."
			#dont know why ts does not accept kmc="python3 ...."
                       	echo $sudoPW| sudo sed -i "s|$previous_kmc|alias kmc=\"ts python3 $package_path/$paralell_file\"|g" "$bash_path/.bashrc"   
			echo $sudoPW| sudo source "$bash_path/.bashrc"
		else
			echo "You do not have previous registration of a kmc alias, only adding..."
			#dont know why ts does not accept kmc="python3 ...."
			new_kmc_entry="alias kmc=\"ts python3 $package_path/$paralell_file\""
			echo $sudoPW | sudo echo $new_kmc_entry >> "$bash_path/.bashrc"
			echo $sudoPW | sudo source "$bash_path/.bashrc"
		fi

		echo
		echo "Ok, Local install with alias finished. It is expected that the running script (paralell.py) is located in:"
		echo "$package_path"
	fi

	
	if [ $bash_ans = "n" ]; then
		mv -f $source_dir/$paralell_file $hm/		
		echo "Ok, Local install without alias finished. I moved paralell.py to this folder in case of you wanting to run paralell jobs."
	fi
fi


if [ $ans = "2" ]; then

        echo "Please insert your username (ex.: padawan, john..)"
        read username
        
        
        bash_path="/home/$username"
	package_path="/home/$username/KMC_source"
	
	sed  -i "s|$param_bib|bib_path = \"$package_path/\"|g" $param_file
	sed  -i "s|$paralell_bib|bib_path = \"$package_path/\"|g" "$source_dir/$paralell_file"	
	
	
	mkdir -p $package_path
	mv -f $source_dir/* $package_path




	if grep -q "alias kmc" "/home/$username/.bashrc"
		then

		echo "I found the following entry in .bashrc with kmc keyword:"
		echo
		previous_kmc=$(grep "alias kmc" "$bash_path/.bashrc")
		echo $previous_kmc
		echo
		echo "Overwriting..."
		#dont know why ts does not accept kmc="python3 ...."
		echo $sudoPW | sudo sed -i "s|$previous_kmc|alias kmc=\"ts python3 $package_path/$paralell_file\"|g" "$bash_path/.bashrc"
		echo $sudoPW | sudo source "$bash_path/.bashrc"
	else
		echo "You do not have previous registration of a kmc alias, only adding..."
		#dont know why ts does not accept kmc="python3 ...."
                new_kmc_entry="alias kmc=\"ts python3 $package_path/$paralell_file\""
		echo $sudoPW | echo $new_kmc_entry >> "$bash_path/.bashrc"
		echo $sudoPW | source "$bash_path/.bashrc"
	fi

	echo
	echo "Ok, Global install with alias finished. It is expected that the running script and source files are located in:"
	echo "$package_path"

fi
