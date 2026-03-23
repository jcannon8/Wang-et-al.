# anEdt7.tcl: Check Nck binding to edt 200 ns MD using NckPRS.dcd references.
# This uses 25 references for PRS-bound Nck from NckPRS.dcd.
# This only tests Nck collision with protein or membrane.
# Output edtPRS10.$mod.dat
# Derived from anEdt5.tcl 
# vmd -dispdev none on Eire.
# Contacts are heavy atoms less than or equal to 1.0 angstroms. 
set dist 1.0;  # Nck collision distance for heavy atoms
cd /t/tcr16/edt
###########################################################
# Load the Nck-PRS references made by makeNckref.tcl
set nckRef [mol new /t/tcr39an/NckPRS.pdb]
mol addfile /t/tcr39an/NckPRS.dcd waitfor all
# Define reference atomselections
set PRSref [atomselect $nckRef "name CA and chain E and resid 180 to 186"]
set nckSH [atomselect $nckRef "chain K"]
set refAll [atomselect $nckRef all]
set numRef [expr [molinfo $nckRef get numframes]-1]
###########################################################
# Load edt 5.4-205 ns MD trajectories for 11 edt models.
set edtModels "84 171 30 115 35 90 13 41 32 153 198"
foreach mod $edtModels {
# Use stripped topology with just 2448 protein atoms that were saved in trajectory.
# made by anEdt1.tcl
set tcr [mol new /t/tcr16/edt.prot.psf]
# The NPT MD started with edt.*.eq2.
for {set i 2} {$i<=53} {incr i} {
  set mdFile /t/tcr16/edt/edt.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
set out [open edtPRS10.$mod.dat w]
###########################################################
# PRS is Pro180-Pro186 in 6JXR chains E and F (see anEdt2.sh).
# Use resid here because edt.prot.psf supplies chain names and can use chain resid.
set PRS [atomselect $tcr "name CA and chain E and resid 180 to 186"]
# notPRS atoms are more than two residues away from PRS
set notPRS [atomselect $tcr "not hydrogen and not (chain E and resid 178 to 188)"] 
# CA of first chain E CT residue Ala157 used to judge membrane location. 
set Mem [atomselect $tcr "name CA and chain E and resid 157"]
set lastframe [molinfo $tcr get numframes]
for {set f 0} {$f<$lastframe} {incr f} {
  $PRS frame $f
  $notPRS frame $f
  $Mem frame $f
  # Start output line with frame#
  # ^This should be omitted because the frame is output below.
  puts -nonewline $out [format "%4d" $f]
  # Next Nck contacts for each nckRef
  for {set ref 0} {$ref<$numRef} {incr ref} {
    $PRSref frame $ref; $nckSH frame $ref; $refAll frame $ref;
    # defaults in case no Nck binding.
    set RMSD 0; set nckMem -200; set nRef -1
    # Move refAll to bind to ePRS
    set M [measure fit $PRSref $PRS]; # Fit the 7 PRS CA atoms
    $refAll move $M
    # Check contacts with Nck SH3.1
    set nckCon [measure contacts $dist $nckSH $notPRS] 
    set conNum [llen [lindex $nckCon 1]]
    if {$conNum==0} {
      set nRef $ref
      set RMSD [measure rmsd $PRSref $PRS]
      # Save Nck to membrane distance
      set loc [measure minmax $nckSH]
      set minZ [lindex $loc 0 2]
      set memZ [$Mem get z]
      set eNckMem [expr $minZ- $memZ]
      break; # Stop searching once the first Nck binding without contact found.
    }
  }
  puts [format "%6d%4d %7.2f %7.2f" $f $nRef $RMSD $eNckMem]
  puts $out [format "%6d%4d %7.2f %7.2f" $f $nRef $RMSD $eNckMem]
  # output <frame><nRef><eNck binding><eNck membrane>
  # Nck collision with membrane if MemNck < -6.0 (check on this).
}
close $out
mol delete $tcr
}
###########################################################
# Analysis for frequency of Nck binding to edt frames.
# This uses the NckPRS.dcd Nck-PRS references.
tclsh<<"eof"
# Disregard the PRS RMSD
set models {84 171 30 115 35 90 13 41 32 153 198} 
foreach mod $models {
puts $mod
set in [open edtPRS10.$mod.dat r]
set inData [read $in]
close $in
set out [open edtPRS11.$mod.dat w]
# Example output of edtPRS10.$mod.dat. Note repeated frame#
#   0     0   0    2.33   21.55
#   1     1   0    2.60   21.52
set nCon 0; set mCon 0
foreach line [split $inData "\n"] {
  if {$line == ""} break
  incr frame
  set tm [lindex $line 0]
  incr tm
  # Get Nck-PRS ref index
  set nRef [lindex $line 2]
  incr r($nRef); # track how many times each Nck-PRS ref used.
  set RMSD [lindex $line 3]
  set dMem [lindex $line 4]
  # Nck can bind if no contacts with other proteins and 
  # no Nck contact with membrane.
  if {$nRef>-1 && $dMem>-6.0 } {incr nCon}
  # If no protein contact 
  if {$nRef>-1 } {incr mCon}
  # In 10 ns blocks, each frame is 10 ps.
  if {[expr $tm % 1000]==0} {
    # output <time (ns)><percent Nck bound><percent no protein contact>
    puts $out [format "%4d %4.3f %4.3f" [expr $tm/100] [expr $nCon/10.0] [expr $mCon/10.0]]
    set nCon 0; set mCon 0
  }
}
close $out
}
# Output how many times each Nck reference bound.
parray r
puts "Total frames=$frame"
eof
###########################################################
# Plot Nck binding frequency for eleven edt models in multiplot
# edt models incluster prevalence order:
models="84 171 30 115 35 90 13 41 32 153 198"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out 'edt.nck8.jpg'
# layout: rows, columns
set multiplot layout 3,4 title "{/:Bold Potential Nck binding to CD3εδ using 25 NckPRS references}"
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 6
set xtics out nomirror 0,50,220
set ytic out nomirror offset 1,0
set xrange [0:210]
set yrange [0:100]
set key bottom center textcolor variable samplen -1 font "arial.ttf,12"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
do for [n in "$models"] {
  set title sprintf("edt.%s",n) offset 0,-1
  if (n==93) {set ylabel "{/:Bold Potential Nck binding frequency}" offset 3.4,0}
  else {set ylabel " "}
#  if (n>56) {set xlabel "{/:Bold Time (ns)}"}
#  else {unset xlabel}
  name2 = sprintf("edtPRS11.%s.dat",n)
  plot name2 using 1 : 2 with linespoints ls 2 t "pro + mem",\
  '' using 1 : 3 with linespoints ls 1 t "pro only"
}
quit
eof
###########################################################
# Here update below
# Summary of Nck binding frequency for whole edt and eds ensembles.
edtModels="84 171 30 115 35 90 13 41 32 153 198"
edsModels="154 16 75 174 147 20 176 197 62 107 201"
rm -f temp1 temp2
for mod in $edtModels; do
echo -n "/t/tcr16/edt/edtPRS11.$mod.dat " >> temp1
echo $mod
done
paste $(cat temp1) > temp2
# 11 edt models
awk '{sum=$2+$5+$8+$11+$14+$17+$20+$23+$26+$29+$32;\
printf("%d %f\n",$1,sum/11)}' temp2>edt.sum2.dat
#
rm -f temp1 temp2
for mod in $edsModels; do
echo -n "/t/tcr16/eds/edsPRS21.$mod.dat " >> temp1
echo $mod
done
paste $(cat temp1) > temp2
# 11 eds models
awk '{sum=$2+$5+$8+$11+$14+$17+$20+$23+$26+$29+$32;\
printf("%d %f\n",$1,sum/11)}' temp2>eds.sum2.dat
rm -f temp1 temp2
###########################################################
# Summary plot of edt and eds ensembles
# This needs to have standard deviation added
gnuplot<<eof
set term jpeg font "arial.ttf,18" 
set out 'edteds.nck8.jpg'
# layout: rows, columns
set title "{/:Bold Nck binding to CD3εδ}"
set border 3
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set xrange [0:210]
set yrange [0:100]
set xlabel "{/:Bold Time (ns)}" offset 0,0.5
set ylabel "{/:Bold Potential Nck binding frequency}" offset 2,0
set key left horizontal textcolor variable samplen -1 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
plot 'edt.sum2.dat' using 1:2 with linespoints ls 1 t "edt",\
'eds.sum2.dat' using (\$1+100):2 with linespoints ls 2 t "eds"
quit
eof
###########################################################
# Use R for statistics
read.table("edt.sum2.dat")->a
mean(a[,2])
read.table("../fgs/fgt.sum2.dat")->b
mean(b[10:21,2])
read.table("../fgs/fgs.sum2.dat")->c
mean(c[10:21,2])

