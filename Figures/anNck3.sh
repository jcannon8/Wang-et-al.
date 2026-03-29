# anNck3.sh: Nck edn interaction (Fig 2A-D).
# Build edr Nck-PRS interaction restraints.
# Analyze Nck-PRS interaction in edn and edm ensembles and edn MD.
# This makes plots for figure and Amber restraint file, edr.RST.
# Using 7-27-24 criteria for edn1E.
cd /t/tcr16/edn
# Use vmd -dispdev none on Eire. 
# Remember VMD residue and index are zero-based. 
# Amber residue and index are one-based.
# This for cd3edm.ref1.pdb model 
set tcr3 [mol new /t/tcr16/charm-gui8/amber/step5_input.psf ]
# This file has chain order e,d,k
mol addfile /t/tcr16/charm-gui8/amber/cd3edm.ref1.pdb waitfor all
# Get atom indices for cpptraj and pmemd restraints.
# Those indices are zero-based, add one to get one-based used in
# cpptraj and Amber restraints.
# 1. PRS Pro656 O – Nck Tyr82 HH
set PRS1 [atomselect $tcr3 "segid PROA and resid 656 and name O"]
set Nck1 [atomselect $tcr3 "segid PROC and resid 82 and name HH"]
$PRS1 get {residue name}; $Nck1 get {residue name}
# 2. PRS Asn660 O – Nck Trp65 HE1
set PRS2 [atomselect $tcr3 "segid PROA and resid 660 and name O"]
set Nck2 [atomselect $tcr3 "segid PROC and resid 65 and name HE1"]
$PRS2 get {residue name}; $Nck2 get {residue name}
# 3. Pro655 CA, CB, CG, CD COM Nck Trp65 CE2, CZ2, CH2, CZ3, CE3, CD2 COM
set PRS3 [atomselect $tcr3 "segid PROA and resid 665 and name CA CB CG CD"]
set Nck3 [atomselect $tcr3 "segid PROC and resid 65 and name CE2 CZ2 CH2 CZ3 CE3 CD2"]
$PRS3 get {residue name}; $Nck3 get {residue name}
# 4. PRS Tyr663 HH – Nck Gln46 OE1
set PRS4 [atomselect $tcr3 "segid PROA and resid 663 and name HH"]
set Nck4 [atomselect $tcr3 "segid PROC and resid 46 and name OE1"]
$PRS4 get {residue name}; $Nck4 get {residue name}
# 5. PRS Tyr663 CG, CD1, CE1, CZ, CE2, CD2 COM – Nck Trp65 CZ2, CH2, CZ3, CE3, CD2, CE2 
# and Tyr77 CG, CD1, CE1, CZ, CE2, CD2 
set PRS5 [atomselect $tcr3 "segid PROA and resid 663 and name CG CD1 CE1 CZ CE2 CD2"]
set Nck5a [atomselect $tcr3 "segid PROC and resid 65 and name CZ2 CH2 CZ3 CE3 CD2 CE2"]
set Nck5b [atomselect $tcr3 "segid PROC and resid 77 and name CG CD1 CE1 CZ CE2 CD2"]
$PRS5 get {residue name}; $Nck5a get {residue name}; $Nck5b get {residue name}
# hydrogen bond between the Val161 carbonyl oxygen and the Asn163 amide proton 
set hb1 [atomselect $tcr3 "segid PROA and resid 658 and name O"]
set hb2 [atomselect $tcr3 "segid PROA and resid 660 and name HN"]
$hb1 get {residue name}; $hb2 get {residue name};
###########################################################
# Use cpptraj to make sure masks work on cd3edm.ref1.pdb.
out=cd3edm.dis1.dat
cpptraj<<eof
parm /t/tcr16/charm-gui8/amber/cd3edm.ref1.pdb
trajin /t/tcr16/charm-gui8/amber/cd3edm.ref1.pdb
# The following works for cd3edm.ref1.pdb with chain edk order
distance PN1 :57@O :209@HH out $out
distance PN2 :61@O :192@HE1 out $out
distance PN3 :66@CD,CA,CB,CG :192@CE2,CD2,CE3,CZ3,CZ2,CH2 out $out
distance PN4 :64@HH :173@OE1 out $out
# This next one has atoms in two residues
distance PN5 :64@CG,CD1,CE1,CZ,CD2,CE2 \
:192@CE2,CD2,CE3,CZ3,CZ2,CH2,:204,CG,CD1,CE1,CZ,CD2,CE2 out $out
distance HB12 :59@O :61@HN out $out
go
eof
#
cat cd3edm.dis1.dat
#Frame            PN1          PN2          PN3          PN4          PN5         HB12
       1       1.7894       2.2350      13.0722       1.9352       5.2351       4.4885
###########################################################
# Get restraints for edn, edm models, which have dek chain order
# These have 6JXR residue numbers.
set edn [mol new /t/tcr16/edn/edn.1.pdb ] 
[atomselect $edn "name CA and chain E and resid 180 to 188"] get {resname resid}
# {PRO 180} {PRO 181} {PRO 182} {VAL 183} {PRO 184} {ASN 185} {PRO 186} {ASP 187} {TYR 188}
# 655        656         657     658      659        660       661       662       663
# 1. PRS Pro181 O – Nck Tyr82 HH
set PRS1 [atomselect $edn "chain E and resid 181 and name O"]
set Nck1 [atomselect $edn "chain K and resid 82 and name HH"]
$PRS1 get {residue name}; $Nck1 get {residue name}
# 2. PRS Asn185 O – Nck Trp65 HE1
set PRS2 [atomselect $edn "chain E and resid 185 and name O"]
set Nck2 [atomselect $edn "chain K and resid 65 and name HE1"]
$PRS2 get {residue name}; $Nck2 get {residue name}
# 3. Pro180 CA, CB, CG, CD COM Nck Trp65 CE2, CZ2, CH2, CZ3, CE3, CD2 COM
set PRS3 [atomselect $edn "chain E and resid 180 and name CA CB CG CD"]
set Nck3 [atomselect $edn "chain K and resid 65 and name CE2 CZ2 CH2 CZ3 CE3 CD2"]
$PRS3 get {residue name}; $Nck3 get {residue name}
# 4. PRS Tyr188 HH – Nck Gln46 OE1
set PRS4 [atomselect $edn "chain E and resid 188 and name HH"]
set Nck4 [atomselect $edn "chain K and resid 46 and name OE1"]
$PRS4 get {residue name}; $Nck4 get {residue name}
# 5. PRS Tyr188 CG, CD1, CE1, CZ, CE2, CD2 COM – Nck Trp65 CZ2, CH2, CZ3, CE3, CD2, CE2 
# and Tyr77 CG, CD1, CE1, CZ, CE2, CD2 
set PRS5 [atomselect $edn "chain E and resid 188 and name CG CD1 CE1 CZ CE2 CD2"]
set Nck5a [atomselect $edn "chain K and resid 65 and name CZ2 CH2 CZ3 CE3 CD2 CE2"]
set Nck5b [atomselect $edn "chain K and resid 77 and name CG CD1 CE1 CZ CE2 CD2"]
$PRS5 get {residue name}; $Nck5a get {residue name}; $Nck5b get {residue name}
# hydrogen bond between the Val161 carbonyl oxygen and the Asn163 amide proton 
set hb1 [atomselect $edn "chain E and resid 183 and name O"]
set hb2 [atomselect $edn "chain E and resid 185 and name HN"]
$hb1 get {residue name}; $hb2 get {residue name};
# PRS Asp187 to Nck Lys37 (called Lys36 in Takeuchi in error. 
set N187 [atomselect $edn "name CG and chain E and resid 187"]
set K36 [atomselect $edn "name NZ and chain K and resid 37"]
$N187 get {residue name}; $K36 get {residue name}; # one-based residues 137, 164
###########################################################
# Use bash
# Use indices to make PN1-5 restraints
cat <<eof> edr.RST
# PN1 maximum 10.0
&rst
    iat=2013, 3321,
    r1=0, r2=0, r3=10.0, r4=12.0
    rk2=10.0, rk3=10.0,
/
# PN2 maximum 12.0
&rst
    iat=2071, 3037,
    r1=0, r2=0, r3=12.0, r4=15.0
    rk2=10.0, rk3=10.0,
/
# PN3 maximum 20.0
&rst
    iat=-1, -1,
    igr1=1987,1990,1992,1995,
    igr2=3038,3039,3040,3042,3044,3046,
    r1=0, r2=0, r3=20.0, r4=21.0
    rk2=10.0, rk3=10.0,
/
/
# PN5 maximum 25.0
&rst
    iat=-1, -1,
    igr1=2105,2106,2108,2110,2113,2115,
    igr2=3038,3039,3040,3042,3044,3046,3238,3239,3241,3243,3246,3248,
    r1=0, r2=0, r3=25.0, r4=27.0
    rk2=10.0, rk3=10.0,
/
eof

###########################################################
# Use cpptraj to make sure masks work on ednEns and edmEns ensemble.
out=edm.dis.dat
cpptraj<<eof
parm edn/edn.1.pdb
# Only protein atoms saved in cdf
parmstrip !@1-3405
# ednEns or edmEns ensemble
trajin edmEns.dcd
distance PN1 :131@O :209@HH out $out
distance PN2 :135@O :192@HE1 out $out
distance PN3 :130@CD,CA,CB,CG :192@CE2,CD2,CE3,CZ3,CZ2,CH2 out $out
distance PN4 :138@HH :173@OE1 out $out
# This next one has atoms in two residues
distance PN5 :138@CG,CD1,CE1,CZ,CD2,CE2 \
:192@CE2,CD2,CE3,CZ3,CZ2,CH2,:204,CG,CD1,CE1,CZ,CD2,CE2 out $out
distance HB12 :133@O :135@HN out $out
go
eof
# Combine Nck-PRS contacts from PRScomp2.tcl with this restraint data.
echo "model con">temp
cat edm.prs2.dat>>temp
paste edm.dis.dat temp>temp2
less temp2
###########################################################
# USe bash, ccptraj
# Analyze restraint distances and Nck-PRS LIE in edn 100 ns MD
# 115 105 126 11 194 157 9 152 181 106 103 56 54  
for mod in 115 105 126 11 194 157 9 152 181 106 103 56 54; do
name=edn.$mod
echo $name
# Build collection of trajectory inputs from eq14 to eq33.
rm -f trajx
for ((j=14; j<=33;j++)); do
echo trajin /t/tcr16/edn/edn.$mod.eq$j.cdf >> trajx
done
# Run ccptraj to get restrain values
out=$name.dist.dat
out2=$name.lie.dat
cpptraj<<eof
# LIE needs a topology with force field parameters, a pdb will not work 
parm /t/tcr16/edn/edn.$mod.top
# Only protein atoms saved in cdf
parmstrip !@1-3405
# Load trajectories like "source"
readinput trajx
distance PN1 :131@O :209@HH out $out
distance PN2 :135@O :192@HE1 out $out
distance PN3 :130@CD,CA,CB,CG :192@CE2,CD2,CE3,CZ3,CZ2,CH2 out $out
distance PN4 :138@HH :173@OE1 out $out
# This next one has atoms in two residues
distance PN5 :138@CG,CD1,CE1,CZ,CD2,CE2 \
:192@CE2,CD2,CE3,CZ3,CZ2,CH2,:204,CG,CD1,CE1,CZ,CD2,CE2 out $out
# PRS Val183 carbonyl to Asn185 amide
distance HB12 :133@O :135@HN out $out
# PRS Asp187 to Nck Lys36
distance salt :137@CG :164@NZ out $out
# Linear interaction energy between Nck and chain E (PRS or all) and chain D
# Nck (chain K) residues 158-213 (one-based)
# CD3e PRS residues 130-138 (one-based), CD3e residues 75-157 (one-based)
# CD3d (chain D) residues 1-74 (one-based)
lie PRS :130-138 :158-213 out $out2
lie CD3e :75-157 :158-213 out $out2
lie CD3d :1-74 :158-213 out $out2
go
eof
done
rm -f trajx
###########################################################
# Figure 2B
# Plot restraint distances for 13 edn 100 ns MD
# See https://makemeengr.com/gnuplot-apply-colornames-from-datafile/ for colors
# Use the hexidecimal colorspec as written in this legend table.
cd /t/tcr16/edn
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edn.dis2.jpg'
set multiplot layout 3,5 title "{/:Bold CD3εδ^N Nck-PRS distances}"
set border 3
set rmargin 1.5
set lmargin 3
set xtics out nomirror 0,50,120
set xrange [0:110]
set ytic out nomirror 0,10,40 offset 1,0
set yrange [0:40]
set key Left reverse textcolor variable samplen -1 font "arial.ttf,18"
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 3 pt 7 pi -1 ps 1.0
#
L=0.62
do for [n in "115 105 126 11 194 157 9 152 181 106 103 56 54" ] {
  set title sprintf("edn.%s",n) offset 0,-1
  filename = sprintf("edn.%s.dist.dat",n)
  plot filename using (\$1/100) : 2 with lines ls 1 title \
  "PRS Pro181 O – Nck Tyr82 HH" at L,0.25,\
  '' using (\$1/100) : 3 with lines ls 2 title \
  "PRS Asn185 O – Nck Trp65 HE1" at L, 0.22,\
  '' using (\$1/100) : 4 with lines ls 3 title \
  "PRS Pro180 - Nck Trp65" at L, 0.19,\
  '' using (\$1/100) : 6 with lines ls 4 title \
  "PRS Tyr188 - Nck Trp65,Tyr77" at L, 0.16,\
  '' using (\$1/100) : 8 with lines ls 5 title \
  "PRS Asp187 - Nck Lys37" at L, 0.13
}
quit
eof
###########################################################
# Figure 2A
# Plot all Nck LIE for 13 edn 
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edn1b.lie.jpg'
set multiplot layout 3,5 title "{/:Bold CD3εδ^N Nck-CD3εδ interaction energy}"
set border 3
set rmargin 1.5
set lmargin 4.2
set xtics out nomirror 0,50,120
set xrange [0:110]
set ytic out -200,50,0 nomirror offset 1,0
set yrange [-200:10]
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
L=0.62
do for [n in "115 105 126 11 194 157 9 152 181 106 103 56 54" ] {
  set title sprintf("edn.%s",n) offset 0,-1
  fname = sprintf("edn.%s.lie.dat",n)
  plot fname u (\$1/100) : (g(\$2,\$3)) with lines ls 2 t "Nck-PRS" at L,0.19,\
  '' u (\$1/100) : (g(\$4,\$5)) with lines ls 1 t "Nck-CD3ε" at L, 0.16,\
  '' u (\$1/100) : (g(\$6,\$7)) with lines ls 3 t "Nck-CD3εδ" at L, 0.13
}
quit
eof
###########################################################
# Check CD3 PRS dihedral angles for Nck binding criteria in edn
# vmd -dispdev none on Eire.
cd /t/tcr16/eds
# Using wrapped angles, these limits prevent Nck binding
# P182 ψ <50°, V183 φ<240°, P184 ψ <50° 
proc wrapPhi {a} {
  return [expr $a<0?$a+360:$a]
}
proc wrapPsi {a} {
  return [expr $a<-120?$a+360:$a]
}
set ednModels "115 105 126 11 194 157 9 152 181 106 103 56 54"
foreach mod $ednModels {
set out [open ednPRS16.$mod.dat w]
# Use stripped topology with just 3405 protein atoms that were saved in trajectory.
set tcr [mol new /t/tcr16/edn.prot.psf]
#[atomselect $tcr "name CA"] get {chain resname resid}
# ^To check resid
# Build collection of trajectory inputs from eq14 to eq33, 100 ns total.
for {set i 14} {$i<=33} {incr i} {
  set mdFile /t/tcr16/edn/edn.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
# Use edr resid from /t/tcr16/edn.prot.psf
set P182 [atomselect $tcr "chain E and name CA and resid 182"]
set V183 [atomselect $tcr "chain E and name CA and resid 183"]
set P184 [atomselect $tcr "chain E and name CA and resid 184"]
set lastframe [molinfo $tcr get numframes]
for {set f 0} {$f<$lastframe} {incr f} {
  set badA 0; set bad182 0; set bad183 0; set bad184 0
  $P182 frame $f; $V183 frame $f; $P184 frame $f
  # These are the 7-27-24 criteria
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
# Below are the 7-18-24 criteria
if{0}
  set P182psi [$P182 get psi] 
  if {[wrapPsi $P182psi] < 50} {incr badA; incr bad182}
  set P183phi [$V183 get phi]
  if {[wrapPhi $P183phi] < 240} {incr badA; incr bad183}
  set P184psi [$P184 get psi]
  if {[wrapPsi $P184psi] < 50} {incr badA; incr bad184}
} 
###########################################################
# Figure 2D
# Plot PRS angle criteria for edn
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edn.ang.jpg'
set multiplot layout 3,5 title "{/:Bold CD3εδ^N PRS angle violations}"
set border 3
set rmargin 1.5
set lmargin 4.3
set xtics out nomirror 0,50,120
set xrange [0:110]
set ytic out 0,1,2 nomirror offset 1,0
set yrange [-0.5:2.5]
set grid ytics
# These all have thin lines.
set style line 1 lc rgb "red" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 1 pt 7 pi -1 ps 1.0
#
do for [n in "115 105 126 11 194 157 9 152 181 106 103 56 54" ] {
  set title sprintf("edn.%s",n) offset 0,-1
  fname = sprintf("ednPRS16.%s.dat",n)
  plot fname u (\$1/100):2 with points ls 2 t ""
}
quit
eof
###########################################################
cat ednPRS16.*.dat > ednPRS16.sum.dat
# use R to get summary of violations
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
read.table("ednPRS16.sum.dat")->a
length(a[,1])
[1] 129557
> sum(a[,2])/length(a[,1])
[1] 0.4098273
> sum(a[,3])/length(a[,1])
[1] 0.1555223
> sum(a[,4])/length(a[,1])
[1] 0.2543051
###########################################################
# LIE distribution for edn ensemble. NOT in final form!
# Concatenate data from all models in ensemble..
cat edn.*.lie.dat> ednEns.lie.dat
# Concatenate only data from models with CD3d <-10.
awk '{if($6<-10 && $7<-10) print $0}' edn.*.lie.dat>ednEns2.lie.dat
#
gnuplot<<eof
#set term jpeg font "arial.ttf,18" 
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edn2c.lie.jpg'
set title "{/:Bold Edn Nck-CD3εδ interaction energy\n}{/:Bold distribution}"
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
r=0.001
plot "ednEns2.lie.dat" u (bin((g(\$2,\$3)),1.0)):(r) s f w lines ls 2 t "Nck-PRS",\
'' u (bin((g(\$4,\$5)),1.0)):(r) s f w lines ls 1 t "Nck-CD3ε",\
"ednEns2.lie.dat" u (bin((g(\$6,\$7)),1.0)):(r) s f w lines ls 3 t "Nck-CD3δ"
quit
eof
###########################################################
# Analyze Nck and PRS RMSD in edn 100 ns MD
# 115 105 126 11 194 157 9 152 181 106 103 56 54  
for mod in 115 105 126 11 194 157 9 152 181 106 103 56 54; do
name=edn.$mod
echo $name
# Build collection of trajectory inputs
rm -f trajx
for ((j=14; j<=33;j++)); do
echo trajin /t/tcr16/edn/edn.$mod.eq$j.cdf >> trajx
done
# Run ccptraj to get restraint distances
out=$name.dist.dat
cpptraj<<eof
# Use the edm topology
parm /t/tcr16/edm/edm.$mod.top
# Only 3405 protein atoms in trajectory
parmstrip !@1-3405
readinput trajx
# Nck reference is frame 9803 of 2jxb3 MD.
parm /t/tcr16/eds/Nck9803.pdb [ref0]
# Use ref0 topology and name reference "ref2".
reference /t/tcr16/eds/Nck9803.pdb parm [ref0] [ref2]
# PRS reference is also frame 9803
parm /t/tcr16/eds/NckPRS9803.pdb [ref3]
reference /t/tcr16/eds/NckPRS9803.pdb parm [ref3] [ref4]
# Note order of masks: first mask is for trajectory, second for reference.
rmsd Nck ref [ref2] ^3@CA @CA out nck.rms.$mod.dat
# PRS Amber residue numbers (one-based) determined via VMD above
rmsd PRS ref [ref4] :130-138@CA :6-14@CA out nck.rms.$mod.dat
eof
done
rm -f trajx
###########################################################
# Figure 2C
# Plot Nck and PRS RMSD for edn models
models="115 105 126 11 194 157 9 152 181 106 103 56 54"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edn3.rms.jpg'
# layout: rows, columns
set multiplot layout 3,5 title "{/:Bold CD3εδ^N Nck and PRS RMSD}"
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 4
set xtics out nomirror 0,50,120
set ytic out nomirror 0,1,4 offset 1,0
set xrange [0:110]
set yrange [0:4]
set key left horizontal textcolor variable samplen -1 font "arial.ttf,12"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
L=0.66
do for [n in "$models"] {
  set title sprintf("edn.%s",n) offset 0,-1
#  if (n==147) {set ylabel "{/:Bold RMSD}" offset 4.2,0}
#  else {set ylabel " "}
#  if (n>56) {set xlabel "{/:Bold Time (ns)}"}
#  else {unset xlabel}
  name2 = sprintf("nck.rms.%s.dat",n)
  plot name2 using (\$1/100):2 with lines ls 1 t "Nck" at L,0.16,\
  '' using (\$1/100):3 with lines ls 2 t "PRS" at L,0.19
}
quit
eof
###########################################################
# Use vmd -dispdev none
# Get restraints for fgn, fgr, and fgm models, which have FGK chain order
# These have 6JXR residue numbers.
set fgm [mol new /t/tcr16/fgm/fgm.1.pdb ] 
[atomselect $fgm "name CA and chain F and resid 180 to 188"] get {resname resid}
# {PRO 180} {PRO 181} {PRO 182} {VAL 183} {PRO 184} {ASN 185} {PRO 186} {ASP 187} {TYR 188}
# 655        656         657     658      659        660       661       662       663
# 1. PRS Pro181 O – Nck Tyr82 HH
set PRS1 [atomselect $fgm "chain F and resid 181 and name O"]
set Nck1 [atomselect $fgm "chain K and resid 82 and name HH"]
$PRS1 get {residue name}; $Nck1 get {residue name}
# 2. PRS Asn185 O – Nck Trp65 HE1
set PRS2 [atomselect $fgm "chain F and resid 185 and name O"]
set Nck2 [atomselect $fgm "chain K and resid 65 and name HE1"]
$PRS2 get {residue name}; $Nck2 get {residue name}
# 3. Pro180 CA, CB, CG, CD COM Nck Trp65 CE2, CZ2, CH2, CZ3, CE3, CD2 COM
set PRS3 [atomselect $fgm "chain F and resid 180 and name CA CB CG CD"]
set Nck3 [atomselect $fgm "chain K and resid 65 and name CE2 CZ2 CH2 CZ3 CE3 CD2"]
$PRS3 get {residue name}; $Nck3 get {residue name}
# 4. PRS Tyr188 HH – Nck Gln46 OE1
set PRS4 [atomselect $fgm "chain F and resid 188 and name HH"]
set Nck4 [atomselect $fgm "chain K and resid 46 and name OE1"]
$PRS4 get {residue name}; $Nck4 get {residue name}
# 5. PRS Tyr188 CG, CD1, CE1, CZ, CE2, CD2 COM – Nck Trp65 CZ2, CH2, CZ3, CE3, CD2, CE2 
# and Tyr77 CG, CD1, CE1, CZ, CE2, CD2 
set PRS5 [atomselect $fgm "chain F and resid 188 and name CG CD1 CE1 CZ CE2 CD2"]
set Nck5a [atomselect $fgm "chain K and resid 65 and name CZ2 CH2 CZ3 CE3 CD2 CE2"]
set Nck5b [atomselect $fgm "chain K and resid 77 and name CG CD1 CE1 CZ CE2 CD2"]
$PRS5 get {residue name}; $Nck5a get {residue name}; $Nck5b get {residue name}
# hydrogen bond between the Val161 carbonyl oxygen and the Asn163 amide proton 
set hb1 [atomselect $fgm "chain F and resid 183 and name O"]
set hb2 [atomselect $fgm "chain F and resid 185 and name HN"]
$hb1 get {residue name}; $hb2 get {residue name};
# Use indices to make PN1-5 restraints
$PRS1 get index; $Nck1 get index
$PRS2 get index; $Nck2 get index
$PRS3 get index; $Nck3 get index
$PRS5 get index; $Nck5a get index; $Nck5b get index
# make restraint file
cat <<eof> fgr.RST
# PN1 maximum 10.0
&rst
    iat=929, 3425,
    r1=0, r2=0, r3=10.0, r4=12.0
    rk2=10.0, rk3=10.0,
/
# PN2 maximum 12.0
&rst
    iat=987, 3141,
    r1=0, r2=0, r3=12.0, r4=15.0
    rk2=10.0, rk3=10.0,
/
# PN3 maximum 20.0
&rst
    iat=-1, -1,
    igr1=903,906,908,911,
    igr2=3142,3143,3144,3146,3148,3150,
    r1=0, r2=0, r3=20.0, r4=21.0
    rk2=10.0, rk3=10.0,
/
# PN5 maximum 25.0
&rst
    iat=-1, -1,
    igr1=1021,1022,1024,1026,1029,1031,
    igr2=3142,3143,3144,4146,3148,3150,3342,3343,3345,3347,3350,3352,
    r1=0, r2=0, r3=25.0, r4=27.0
    rk2=10.0, rk3=10.0,
/
eof
###########################################################
