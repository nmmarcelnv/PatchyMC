#!/usr/bin/bash


for T in 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74; do 
	dir=T${T}
	mkdir ${dir}	
	sed 's/TEMP/'${T}'/' <moab.sub >${dir}/moab${T}
	cd ${dir}
	bash moab${T}
	cd ..
done 

