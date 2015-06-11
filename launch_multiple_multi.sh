#!/bin/bash

occ_dir="$PWD"

for species in "Aglais_io" "Aglais_urticae" "Boloria_selene" "Celastrina_argiolus" "Colias_crocea" "Limenitis_camilla" "Lycaena_phlaeas" "Ochlodes_sylvanus" "Papilio_machaon" "Pararge_aegeria" "Pieris_napi" "Pieris_rapae" "Plebejus_argus" "Polygonia_c-album" "Polyommatus_icarus";
 
do

	cp $occ_dir/*.R $occ_dir/$species
	cp $occ_dir/sub.sh $occ_dir/$species
	cp $occ_dir/reto.sh $occ_dir/$species
	cd $occ_dir/$species
	sh sub.sh
	cd ..

done

