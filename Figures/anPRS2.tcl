# anPRS2.tcl: Analyze PRS dihedral angles in 2jxb and 2k4f (Fig 4 A,B,E,F).
cd /t/tcr16/eds
# Get backbone dihedral angles with VMD.
# vmd -dispdev none
# Wrap angles as decribed by Hollingsworth and Karplus 2010
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3061398/
# Wrapped Ramachandran plots reset phi<0 to phi+360 and psi< -120 to psi+360.
# Then plot -120 to 240 psi on y-axis, 0 to 360 phi on x-axis.
proc phiOut {phi outFile last} {
# Output one line of phi angles with wrapping.
  for {set j 0} {$j<=8} {incr j} {
    set a [lindex $phi $j]
    set b [expr $a<0?$a+360:$a]
    puts -nonewline $outFile [format "%6.3f " $b]
  }
  puts $outFile $last
}
proc psiOut {psi outFile last} {
# Output one line of psi angles with wrapping.
  for {set j 0} {$j<=8} { incr j} {
    set a [lindex $psi $j]
    set b [expr $a<-120?$a+360:$a]
    puts -nonewline $outFile [format "%6.3f " $b]
  }
  puts $outFile $last
}
# 2jxb.[phi,psi].dat get angle data from 2jxb.pdb, 2jxb3.c0, and 5qu2.pdb  
set 2jxb [mol new /home/jcannon/tcr2/cd3e/2jxb.pdb waitfor all]
set 2jxb3 [mol new /home/jcannon/tcr2/cd3e/2jxb3.c0.pdb waitfor all]
set nFrames [molinfo $2jxb get numframes]
[atomselect $2jxb "name CA"] get {resname resid}
# PRS is resid 6-14 in 2jxb
set 5qu2 [mol new /t/tcr16/eds/5qu2.pdb waitfor all]
[atomselect $5qu2 "name CA"] get {chain resname resid}
# There are two Nck-PRS complexes in the 5qu2 asymmetric unit, chains D and E.
# PRS is resid 180-188 in 5qu2 
set out1 [open "2jxb.phi.dat" w]
set out2 [open "2jxb.psi.dat" w]
for {set i 0} {$i<=$nFrames} {incr i} {
  set phi [[atomselect $2jxb "name CA and resid 6 to 14" frame $i] get phi]
  phiOut $phi $out1 1
  set psi [[atomselect $2jxb "name CA and resid 6 to 14" frame $i] get psi]
  psiOut $psi $out2 1
}
# Add 5qu2 chain D angles
set psi [[atomselect $5qu2 "name CA and chain D and resid 180 to 188"] get phi]
phiOut $phi $out1 2
set psi [[atomselect $5qu2 "name CA and chain D and resid 180 to 188"] get psi]
psiOut $psi $out2 2
# Add 5qu2 chain E angles
set phi [[atomselect $5qu2 "name CA and chain E and resid 180 to 188"] get phi]
phiOut $phi $out1 2
set psi [[atomselect $5qu2 "name CA and chain D and resid 180 to 188"] get psi]
psiOut $psi $out2 2
# Add 2jxb3.c0 angles
set phi [[atomselect $2jxb3 "name CA and resid 6 to 14"] get phi]
phiOut $phi $out1 3
set psi [[atomselect $2jxb3 "name CA and resid 6 to 14"] get psi]
psiOut $psi $out2 3
#
close $out1
close $out2
#############################
set 2k4f [mol new /home/jcannon/DYFZG02/bifab/2k4f.pdb waitfor all]
set nFrames [molinfo $2k4f get numframes]
[atomselect $2k4f "name CA"] get {resname resid}
[atomselect $2k4f "name CA and resid 30 to 38"] get {resname resid}
# PRS is resid 30-38 in 2k4f
set out1 [open "2k4f.phi.dat" w]
set out2 [open "2k4f.psi.dat" w]
for {set i 0} {$i<=$nFrames} {incr i} {
  set phi [[atomselect $2k4f "name CA and resid 30 to 38" frame $i] get phi]
  phiOut $phi $out1 1
  set psi [[atomselect $2k4f "name CA and resid 30 to 38" frame $i] get psi]
  psiOut $psi $out2 1
}
close $out1
close $out2
###########################################################
# Plot 2jxb and 2k4f angles
# Fig 4 A,B,E,F four panels
gnuplot<<"eof"
#set term jpeg font "arial.ttf,18" size 1280,960
set term jpeg font "arial.ttf,18" size 1000,750
set out 'pdb2.dih.jpg'
# layout: rows, columns
set multiplot layout 2,2 
set border 3
set xtic nomirror rotate by -45 offset -1
set ytic out nomirror
# linetype controls the color of points rather than line style
set linetype 1 lc rgb "blue" 
set linetype 2 lc rgb "red" 
set linetype 3 lc rgb "green" 
set xrange [0:9.5]
set grid ytics
# offset points other than model 1
jit(x,m)=m==1?x:x+0.1
# Do not label offset angle name
lab(s,m)=m==1?s:""
# plot-specific parameters for 2jxb psi
root="2jxb"
set ylabel "{/:Bold ψ angle (degrees)}" offset 1,0
set title "{/:Bold Nck-bound in 2JXB and 5QU2}"
#set label 1 "{/:Bold Nck-bound CD3ε PRS in 2JXB and 5QU2}" at graph 0,1
ext="psi"
set yrange [-120:240]
fname=sprintf("%s.%s.dat",root,ext)
# Each column has a different angle, each row has a different model.
# Tenth colum has 1-3 for point color
plot fname using (jit(1,$10)):1:10:xtic(lab("P180",$10)) with points lc var pt 7 t "",\
'' u (jit(2,$10)):2:10:xtic(lab("P181",$10)) with points lc var pt 7 t "",\
'' u (jit(3,$10)):3:10:xtic(lab("P182",$10)) with points lc var pt 7 t "",\
'' u (jit(4,$10)):4:10:xtic(lab("V183",$10)) with points lc var pt 7 t "",\
'' u (jit(5,$10)):5:10:xtic(lab("P184",$10)) with points lc var pt 7 t "",\
'' u (jit(6,$10)):6:10:xtic(lab("N185",$10)) with points lc var pt 7 t "",\
'' u (jit(7,$10)):7:10:xtic(lab("P186",$10)) with points lc var pt 7 t "",\
'' u (jit(8,$10)):8:10:xtic(lab("D187",$10)) with points lc var pt 7 t "",\
'' u (jit(9,$10)):9:10:xtic(lab("Y188",$10)) with points lc var pt 7 t ""
# Second plot top right for 2jxb phi
set ylabel "{/:Bold φ angle (degrees)}" offset 1,0
#set title "{/:Bold 2jxb φ dihedral angles}"
unset label 1
ext="phi"
set yrange [0:360]
fname=sprintf("%s.%s.dat",root,ext)
# repeat the last plot command with new arguments
replot
# Third plot bottom left for 2k4f psi
root="2k4f"
set ylabel "{/:Bold ψ angle (degrees)}" offset 1,0
set title "{/:Bold Nck-free in 2K4F}"
#set label 1 "{/:Bold Nck-free CD3ε PRS in 2K4F}" at graph 0,1
ext="psi"
set yrange [-120:240]
fname=sprintf("%s.%s.dat",root,ext)
replot
# Fourth plot bottom right for 2k4f phi
set ylabel "{/:Bold φ angle (degrees)}" offset 1,0
#set title "{/:Bold 2k4f φ dihedral angles}"
unset label 1
ext="phi"
set yrange [0:360]
fname=sprintf("%s.%s.dat",root,ext)
replot
#
quit
eof

