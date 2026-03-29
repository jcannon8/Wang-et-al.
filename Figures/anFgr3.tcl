# anFgr3.tcl: Analyze Nck-PRS interactions in fgr (S4 fig, Fig 3B).
# This makes the four plots for (S4 fig.
cd /t/tcr16/fgs
# vmd -dispdev none on Eire.
# Get residue numbers for PRS in fgn model
set tcr [mol new /t/tcr16/fgn.prot.pdb waitfor all]
[atomselect $tcr "chain F and name CA"] get {resname resid}
# PRS is Pro180-Tyr188
[atomselect $tcr "chain F and name CA and resid 180 to 188"] get {resname resid residue}
# ^That shows PRS are residues 57-65 (zero-based), 58-66 (one-based).
[atomselect $tcr "chain F and name CA"] get {resname resid residue}
# ^That shows CD3 chain F are residues 0-207 (zero-based), 1-85 (one-based).
[atomselect $tcr "chain K and name CA"] get {resname resid residue}
# ^That shows Nck are residues 160-215 (zero-based), 161-216 (one-based).
[atomselect $tcr "chain G and name CA"] get {resname resid residue}
# ^That shows CD3 chain G are residues 85-159 (zero-based), 86-160 (one-based).
###########################################################
# Analyze Nck-PRS distances in fgr ensemble
# 2 105 177 199 93 10 83 168 142 30 55 154 117 86   
for mod in 2 105 177 199 93 10 83 168 142 30 55 154 117 86 201; do
name=fgr.$mod
echo $name
# Build collection of trajectory inputs
rm -f trajx
# 420 ns
for ((j=13; j<=95;j++)); do
if [ -e /t/tcr16/fgr/fgr.$mod.eq$j.cdf ]; then
echo trajin /t/tcr16/fgr/fgr.$mod.eq$j.cdf >> trajx
fi
done
# Run ccptraj to get restraint values
out=$name.dist.dat
out2=$name.lie.dat
cpptraj<<eof
# The fgr models have the same toplogy as the same fgn model.
parm /t/tcr16/fgn/fgn.$mod.top
# Only protein atoms saved in cdf
parmstrip !@1-3509
# Load trajectories like "source"
readinput trajx
distance PN1 :59@O :212@HH out $out
distance PN2 :63@O :195@HE1 out $out
distance PN3 :58@CD,CA,CB,CG :195@CE2,CD2,CE3,CZ3,CZ2,CH2 out $out
distance PN4 :66@HH :176@OE1 out $out
# This next one has atoms in two residues
distance PN5 :66@CG,CD1,CE1,CZ,CD2,CE2 \
:195@CE2,CD2,CE3,CZ3,CZ2,CH2,:207,CG,CD1,CE1,CZ,CD2,CE2 out $out
distance HB12 :61@O :63@HN out $out
# PRS Asp187 to Nck Lys36
distance salt :65@CG :167@NZ out $out
# Linear interaction energy between chain F (PRS or all) and Nck
# chain F PRS residues 58-66 (one-based) or chain F residues 1-85 (one-based)
# and Nck (chain K) residues 161-216 (one-based).
lie PRS :58-66 :161-216 out $out2
lie CD3f :1-85 :161-216 out $out2
lie CD3g :86-160 :161-216 out $out2
go
eof
done
rm -f trajx
###########################################################
# Figure S4 Fig panel B
# Plot fgr Nck-PRS distances
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'fgr.dis.jpg'
set multiplot layout 4,4 title "{/:Bold CD3εγ^R Nck-PRS distances}"
set border 3
set rmargin 1.5
set lmargin 3
set tmargin 1
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out nomirror 0,10,40 offset 1,0
set yrange [0:40]
set key Left horizontal textcolor variable samplen -1 font "arial.ttf,14"
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 3 pt 7 pi -1 ps 1.0
L=0.72
do for [n in "2 105 177 199 93 10 83 168 142 30 55 154 117 86 201" ] {
  set title sprintf("fgr.%s",n) offset 0,-1.4
  filename = sprintf("fgr.%s.dist.dat",n)
  plot filename using (\$1/100) : 2 with lines ls 1 title \
  "PRS Pro181 O – Nck Tyr82 HH" at L,0.25,\
  '' using (\$1/100) : 3 with lines ls 2 title \
  "PRS Asn185 O–Nck Trp65 HE1" at L, 0.22,\
  '' using (\$1/100) : 4 with lines ls 3 title \
  "PRS Pro180 - Nck Trp65" at L, 0.19,\
  '' using (\$1/100) : 6 with lines ls 4 title \
  "PRS Tyr188 - Nck Trp65,Tyr77" at L, 0.16,\
  '' using (\$1/100) : 8 with lines ls 5 title \
  "PRS Asp187 - Nck Lys37" at L, 0.13
}
eof
# The "Asn185 O–Nck" needed to be squeezed.
###########################################################
# Figure S4 Fig panel A
# Plot Nck-PRS LIE for fgr
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'fgr.lie.jpg'
set multiplot layout 4,4 title "{/:Bold CD3εγ^R Nck-CD3εγ interaction energy}"
set border 3
set rmargin 1.5
set lmargin 4.3
set tmargin 1
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out -250,50,0 nomirror offset 1,0
set yrange [-250:10]
set grid ytics
set key Left reverse textcolor variable samplen -1 font "arial.ttf,18"
# These all have thin lines.
set style line 1 lc rgb "red" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 1 pt 7 pi -1 ps 1.0
# Function for linear interaction energy (alpha=0.18, beta=0.33)
# Aqvist et al. 1994; Carlsson et al. 2006
g(e,v) = 0.18 * v + 0.33 * e
#
L=0.77
do for [n in "2 105 177 199 93 10 83 168 142 30 55 154 117 86 201" ] {
  set title sprintf("fgr.%s",n) offset 0,-1.4
  fname = sprintf("fgr.%s.lie.dat",n)
  plot fname u (\$1/100) : (g(\$2,\$3)) with lines ls 2 t "Nck-PRS" at L,0.15,\
  '' u (\$1/100) : (g(\$4,\$5)) with lines ls 1 t "Nck-CD3ε" at L, 0.12,\
  '' u (\$1/100) : (g(\$6,\$7)) with lines ls 3 t "Nck-CD3εγ" at L, 0.09
}
quit
eof
###########################################################
# Analyze Nck and PRS RMSD for fgr models.
cd /t/tcr16/fgs
# fgr reference models:
for mod in 2 105 177 199 93 10 83 168 142 30 55 154 117 86 201 ; do
rm -f trajx
#for ((j=13; j<=53;j++)); do
# 420 ns
for ((j=13; j<=95;j++)); do
if [ -e /t/tcr16/fgr/fgr.$mod.eq$j.cdf ]; then
echo trajin /t/tcr16/fgr/fgr.$mod.eq$j.cdf >> trajx
fi
done
cpptraj<<eof
# Use the fgm topology
parm /t/tcr16/fgm/fgm.$mod.top
# Only protein atoms in trajectory
parmstrip !@1-3509
#
readinput trajx
# Reference is frame 9803 of 2jxb3 MD.
parm /t/tcr16/eds/Nck9803.pdb [ref0]
# Use ref0 topology and name reference "ref2".
reference /t/tcr16/eds/Nck9803.pdb parm [ref0] [ref2]
# PRS reference is also frame 9803
parm /t/tcr16/eds/NckPRS9803.pdb [ref3]
reference /t/tcr16/eds/NckPRS9803.pdb parm [ref3] [ref4]
# Note order of masks: first mask is for trajectory, second for reference.
rmsd Nck ref [ref2] ^3@CA @CA out nck.rms.$mod.dat
# PRS Amber residue numbers (one-based) determined via VMD above
rmsd PRS ref [ref4] :58-66@CA :6-14@CA out nck.rms.$mod.dat
eof
done
rm -f trajx
###########################################################
# S4 Fig panel C
# Plot fgr RMSD
models="2 105 177 199 93 10 83 168 142 30 55 154 117 86 201"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'fgr3.rms.jpg'
# layout: rows, columns
set multiplot layout 4,4 title "{/:Bold CD3εγ^R Nck and PRS RMSD}"
set border 3
# margin units are character heights or widths
set rmargin 1.5
set lmargin 4.3
set tmargin 1
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out nomirror 0,1,4 offset 1,0
set yrange [0:4.5]
set key left horizontal textcolor variable samplen -1 font "arial.ttf,12"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
L=0.77
do for [n in "$models"] {
  set title sprintf("fgr.%s",n) offset 0,-1.4
  fname = sprintf("nck.rms.%s.dat",n)
  plot fname using (\$1/100):2 with lines ls 1 t "Nck" at L,0.16,\
  '' using (\$1/100):3 with lines ls 2 t "PRS" at L,0.19
}
quit
eof
###########################################################
# Check CD3 PRS dihedral angles for Nck binding criteria in fgr
# Using wrapped angles (Hollingsworth and Karplus 2010).
proc wrapPhi {a} {
  return [expr $a<0?$a+360:$a]
}
proc wrapPsi {a} {
  return [expr $a<-120?$a+360:$a]
}
# Load fgr 10-418 ns MD trajectories for 15 fgr models (including dissociating fgr.2).
set fgrModels "2 105 177 199 93 10 83 168 142 30 55 154 117 86 201"
foreach mod $fgrModels {
set out [open fgrPRS16.$mod.dat w]
# Use stripped topology with just 3509 protein atoms that were saved in trajectory.
set tcr [mol new /t/tcr16/fgn.prot.psf]
#[atomselect $tcr "name CA"] get {chain resname resid}
# ^To check resid
# The NPT MD started with fgr.*.eq13, 418 ns total.
for {set i 13} {$i<=100} {incr i} {
  set mdFile /t/tcr16/fgr/fgr.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
# Use fgr resid from fgr.prot.psf
set P182 [atomselect $tcr "chain F and name CA and resid 182"]
set V183 [atomselect $tcr "chain F and name CA and resid 183"]
set P184 [atomselect $tcr "chain F and name CA and resid 184"]
set lastframe [molinfo $tcr get numframes]
for {set f 0} {$f<$lastframe} {incr f} {
  set badA 0; set bad182 0; set bad183 0; set bad184 0
  $P182 frame $f; $V183 frame $f; $P184 frame $f
  $P182 frame $f; $V183 frame $f; $P184 frame $f
  # Use the 7-27-24 criteria
  set P182psi [$P182 get psi] 
  if {[wrapPsi $P182psi] < 100} {incr badA; incr bad182}
  set P184psi [$P184 get psi]
  if {[wrapPsi $P184psi] < 100} {incr badA; incr bad184}
  # Output: <frame><total bad angles><bad 182><bad184>
  puts $out [format "%6d%4d%4d%4d" $f $badA $bad182 $bad184]
  # Could output raw or wrapped angles as well.
}
close $out
mol delete $tcr
}
# 7-18-24 angle criteria
if {0} {
  set P182psi [$P182 get psi] 
  if {[wrapPsi $P182psi] < 50} {incr badA; incr bad182}
  set P183phi [$V183 get phi]
  if {[wrapPhi $P183phi] < 240} {incr badA; incr bad183}
  set P184psi [$P184 get psi]
  if {[wrapPsi $P184psi] < 50} {incr badA; incr bad184}
}
###########################################################
# Plot PRS angle violation numbers for fgr
# S4 Fig panel D
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'fgr.ang.jpg'
set multiplot layout 4,4 title "{/:Bold CD3εγ^R PRS angle violations}"
set border 3
set rmargin 1.5
set lmargin 4.3
set tmargin 1
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out 0,1,2 nomirror offset 1,0
set yrange [-0.5:2.5]
set grid ytics
# These all have thin lines.
set style line 1 lc rgb "red" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 1 pt 7 pi -1 ps 1.0
#
do for [n in "2 105 177 199 93 10 83 168 142 30 55 154 117 86 201" ] {
  set title sprintf("fgr.%s",n) offset 0,-1
  fname = sprintf("fgrPRS16.%s.dat",n)
  plot fname u (\$1/100):2 with points ls 2 t ""
}
quit
eof
###########################################################
cat fgrPRS16.*.dat > fgrPRS16.sum.dat
# use R to get summary of violations
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
read.table("fgrPRS16.sum.dat")->a
> sum(a[,2])/length(a[,1])
[1] 0.1801057
> sum(a[,3])/length(a[,1])
[1] 0.0998558
> sum(a[,4])/length(a[,1])
[1] 0.08024988
###########################################################
# LIE distribution for fgr ensemble
# Concatenate data from all models in ensemble.
cat fgr.*.lie.dat> fgrEns.lie.dat
# Concatenate only data from models with CD3γ energies <-10.
awk '{if($6<-10 && $7<-10) print $0}' fgr.*.lie.dat>fgrEns2.lie.dat
# Fig 3B
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 500,375 
set out 'fgr2c.lie.jpg'
#set title "{/:Bold CD3εγ^R Nck-CD3εγ interaction energy\n}{/:Bold distribution}"
set border 3
set xtics out nomirror
set ytics out nomirror offset 1,0 
#set yrange [0:55]
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
# Use linestyle here
set linestyle 1 lc rgb "red" lw 2 
set linestyle 2 lc rgb "blue" lw 2
set linestyle 3 lc rgb "dark-green" lw 2
set style fill transparent pattern 0 noborder
set xlabel "{/:Bold Interaction energy (kcal/mol)}"
set ylabel "{/:Bold Relative frequency}" offset 2,0
# Function for binning
bin(x,s) = s*int(x/s)
# Function for linear interaction energy (alpha=0.18, beta=0.33)
# Aqvist et al. 1994; Carlsson et al. 2006
g(e,v) = 0.18 * v + 0.33 * e
# Plots with smooth, frequency, with lines
# y-axis is relative frequency 
# with r=0.001, max y-axis is 4
# with r=0.004, max y-axis is 30
r=0.004
plot "fgrEns2.lie.dat" u (bin((g(\$2,\$3)),1.0)):(r) s f w lines ls 2 t "Nck-PRS",\
'' u (bin((g(\$4,\$5)),1.0)):(r) s f w lines ls 1 t "Nck-CD3ε",\
"fgrEns2.lie.dat" u (bin((g(\$6,\$7)),1.0)):(r) s f w lines ls 3 t "Nck-CD3γ"
quit
eof

