#!/bin/bash

occ_dir="$PWD"

for species in "Anthocharis_cardamines" "Apatura_iris" "Aphantopus_hyperantus" "Araschnia_levana" "Argynnis_aglaja" "Argynnis_paphia" "Callophrys_rubi" "Carterocephalus_palaemon" "Coenonympha_tullia" "Favonius_quercus" "Gonepteryx_rhamni" "Hesperia_comma" "Hipparchia_semele" "Maniola_jurtina" "Nymphalis_polychloros" "Polyommatus_coridon" "Pyronia_tithonus" "Thymelicus_lineola" "Thymelicus_sylvestris";
 
do

	cp $occ_dir/*.R $occ_dir/$species
	cp $occ_dir/sub.sh $occ_dir/$species
	cp $occ_dir/reto.sh $occ_dir/$species
	cd $occ_dir/$species
	sh sub.sh
	cd ..

done

