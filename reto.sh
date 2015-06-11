#!/bin/sh

cd $TMPDIR
lieu1=$2	

cp -r  $lieu1/* .				
ln -s /afs/in2p3.fr/home/throng/mnhn/bin/R ./R	
cat *.R | R --slave --args $1			

mkdir $lieu1/results
mv jagsoutput.Rdata $lieu1/results
mv *.csv $lieu1/results
