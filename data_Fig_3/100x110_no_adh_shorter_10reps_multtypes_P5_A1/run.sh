#!/bin/bash
if g++ colonies_univ.cpp dynamics_rods_3d.cpp main_nographics.cpp -O3 -w -o $1.exe ; then
#	mkdir $1
#	cp params.h $1/params.h
#	cp run.sh $1/run.sh
#	Eb="1e5"
#	for eat in 10 15 20 25 30 35 40
#	do
#        for death in 20 25 30 35 40
	for p in 5 # 10 20 50
	do
	for a in 0 1 2 5
        do
		VAR=$1"_P"$p"_A"$a
		echo $VAR 
		mkdir $VAR
		cp params.h $VAR/params.h
		cp run.sh $VAR/run.sh
                qsub -N BW$1 submit.sh $1.exe name $VAR Nsamp 10 Stopmode 8 72 P $p A $a
        done
	done
else
echo "error when compiling the program"
fi

