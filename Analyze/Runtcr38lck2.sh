#!/bin/bash
# Runtcr38lck2.sh: Run the tcr38lck2.tcl script for tcr38.*.** ensemble.
rm -f tcr38lck2.dat
for ((j=0;j<=11;j++)); do
for k in 00 01 02 10 11 12 20 21 22; do
name=tcr38.$j.$k 
if [ ! -e /t/tcr38b/$name.psf ]; then
continue
fi 
export mod=$j.$k
vmd -dispdev none -e tcr38lck2.tcl
# Concatentate this model data onto tcr38lck2.dat
cat temp >> tcr38lck2.dat 
#
done
done
