# anTcr38e.tcl: Plot Nck binding to tcr38.* models (S14A fig).
# Initially, this used only protein and membrane collision.
# Now it also includes PRS angle criteria.
# using data from anTcr38d.tcl
cd /t/tcr38an
# Analysis for frequency of Nck binding
# This monitors Nck binding to chain E and F.
tclsh<<"eof"
foreach  i  {0 2 3 5 7 8 9 10 11} {
set mod tcr38.$i
puts $mod
set in [open $mod.PRS2.dat r]
set inData [read $in]
close $in
set out [open $mod.PRS3.dat w]
# Example output of tcr38.0.PRS2.dat
# anTcr38d.tcl output <frame><eNck binding><fNck binding><eNck membrane><fNck membrane>
# Nck collision with membrane if < -6 for eMem -13 for fMem.
# Example of anTcr38d.$mod.$i.dat
#      0  24   0   12.48   9.44
set nConE 0; set nConF 0
foreach line [split $inData "\n"] {
  if {$line == ""} break
  set tm [lindex $line 0]
  if {$tm == 0} continue
  set nRefE [lindex $line 1]
  set nRefF [lindex $line 2]
  set eNckMem [lindex $line 3]
  set fNckMem [lindex $line 4]
  # Nck can bind if no contacts with other proteins and 
  # no Nck contact with membrane.
  if {$nRefE>-1 && $eNckMem>-6.0 } {incr nConE}
  if {$nRefF>-1 && $fNckMem>-13.0 } {incr nConF}
  # In 10 ns blocks, each frame is 10 ps.
  if {[expr $tm % 1000]==0} {
    # output <time (ns)><percent Nck bound to chain E><percent Nck bound to chain F>
    puts $out [format "%4d %4.3f %4.3f" [expr $tm/100] [expr $nConE/10.0] [expr $nConF/10.0]]
    set nConE 0; set nConF 0
  }
}
close $out
}
eof
###########################################################
# Plot Nck binding frequency for Nck binding to chain E and F in tcr38.*.
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out 'tcr38.nck1.jpg'
# layout: rows, columns
set multiplot layout 3,3 
set border 3
set xtics out nomirror 0,50,200
set ytic out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set key center left textcolor variable samplen -1 font "arial.ttf,12"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
do for [n in "0 2 3 5 7 8 9 10 11"] {
  t = sprintf("tcr38.%s",n)
  set title t
  if(n==5) {set ylabel "{/:Bold Potential Nck binding}" offset 2,0}
  else {set ylabel " "} 
  if(n>8) {set xlabel "{/:Bold Time (ns)}"}
  filename = sprintf("tcr38.%s.PRS3.dat",n)
  plot filename using 1:2 with linespoints ls 1 t "CD3εδ (ED)",\
  '' using 1:3 with linespoints ls 2 t "CD3εγ (FG)"
}
quit
eof
###########################################################
# Averages of Nck binding frequency for tcr38 ensemble.
rm -f temp1 temp2
for mod in 0 2 3 5 7 8 9 10 11; do 
echo -n "tcr38.$mod.PRS3.dat " >> temp1
echo $mod
done
paste $(cat temp1) > temp2
# 9 tc38 models
awk '{sum=$2+$5+$8+$11+$14+$17+$20+$23+$26;\
printf("%d %f\n",$1,sum/9)}' temp2>tcr38.sumE.dat
awk '{sum=$3+$6+$9+$12+$15+$18+$21+$24+$27;\
printf("%d %f\n",$1,sum/9)}' temp2>tcr38.sumF.dat
rm -f temp1 temp2
###########################################################
# Summary plot
gnuplot<<eof
set term jpeg font "arial.ttf,18" 
set out 'tcr38.nck2.jpg'
set title "{/:Bold Nck binding to Tcr38}"
set border 3
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set xlabel "{/:Bold Time (ns)}" offset 0,0.5
set ylabel "{/:Bold Potential Nck binding}" offset 2,0
#set key left horizontal textcolor variable samplen -1 
set key left textcolor variable samplen -1
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
plot 'tcr38.sumE.dat' using 1:2 with linespoints ls 1 t "CD3εδ (ED)",\
'tcr38.sumF.dat' using 1:2 with linespoints ls 2 t "CD3εγ (FG)"
quit
eof
###########################################################
###########################################################
# Get PRS angles
# Use the PRS angle criteria
# Modified from anTcr39e.tcl
# vmd -dispdev none
cd /t/tcr38an
# Using wrapped angles (Hollingsworth and Karplus 2010).
proc wrapPhi {a} {
  return [expr $a<0?$a+360:$a]
}
proc wrapPsi {a} {
  return [expr $a<-120?$a+360:$a]
}
###########################################################
# Find a model and load trajectory
#for {set j 0} {$j<=8} {incr j} {
foreach  j {0 2 3 5 7 8 9 10 11} {
set name tcr38.$j
puts $name
# Separate output file for each model
set out [open $name.PRS7.dat w]
# Use protein-only 7668 atom topology
set tcr [mol new /t/tcr38/tcr38.prot.psf]
# The tcr38: 0-211 ns MD eq4 up to 60.
for {set i 4} {$i<=60} {incr i} {
  set mdFile /t/tcr38/$name.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
###########################################################
# Look at critical angles P182psi, P183phi, and P184psi in PRS.
set P182e [atomselect $tcr "chain E and name CA and resid 182"]
set V183e [atomselect $tcr "chain E and name CA and resid 183"]
set P184e [atomselect $tcr "chain E and name CA and resid 184"]
set P182f [atomselect $tcr "chain F and name CA and resid 182"]
set V183f [atomselect $tcr "chain F and name CA and resid 183"]
set P184f [atomselect $tcr "chain F and name CA and resid 184"]
#
set lastframe [molinfo $tcr get numframes]
for {set f 0} {$f<$lastframe} {incr f} {
$P182e frame $f;$V183e frame $f; $P184e frame $f
$P182f frame $f;$V183f frame $f; $P184f frame $f
set eDih 0; set fDih 0
set bad182e 0; set bad183e 0;set bad184e 0;
set bad182f 0; set bad183f 0;set bad184f 0;
# Output wrapped angles
set eP182psi [wrapPsi [$P182e get psi]] 
set eP183phi [wrapPhi [$V183e get phi]]
set eP184psi [wrapPsi [$P184e get psi]]
set fP182psi [wrapPsi [$P182f get psi]] 
set fP183phi [wrapPhi [$V183f get phi]]
set fP184psi [wrapPsi [$P184f get psi]]
# Output of tcr38.*.PRS7.dat: <time>
# ! <eP182psi><eP183phi><eP184psi> | <fP182psi><fP183phi><fP184psi>
puts $out [format "%7d %7.2f %7.2f %7.2f | %7.2f %7.2f %7.2f"\
$f $eP182psi $eP183phi $eP184psi $fP182psi $fP183phi $fP184psi]
if {[expr $f % 1000]==0} {puts $f}
}
close $out
mol delete $tcr
}
###########################################################
# Combine collision data in tcr38.*.PRS2.dat from anTcr38g.tcl
# with angle data tcr38.*.PRS7.dat from above to output tcr38.*.PRS5.dat.
for j in 0 2 3 5 7 8 9 10 11; do
name=tcr38.$j
if [ -e $name.PRS7.dat ]; then
echo $name.PRS2.dat
paste $name.PRS2.dat $name.PRS7.dat > temp
# Collect collision data first, then angle data
awk '{printf("%6d%4d%4d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n",\
$1,$2,$3,$4,$5,$7,$8,$9,$11,$12,$13)}' \
temp>$name.PRS8.dat
fi
# Now $name.PRS8.dat has 
# <frame><eNck protein collision><fNck protein collision>
# <eNck mem dist><fNck mem dist>
# <eP182psi><eP183phi><eP184psi> | <fP182psi><fP183phi><fP184psi>
done
rm -f temp
###########################################################
# Analysis for frequency of Nck binding in 10 ns blocks.
# Use the 7-27-24 PRS angle criteria
# This reports individual criteria as well as cumulative criteria for both Nck. 
tclsh<<"eof"
foreach  j {0 2 3 5 7 8 9 10 11} {
set name tcr38.$j
if {[file exist $name.PRS8.dat]==0} continue
puts $name
set in [open $name.PRS8.dat r]
set inData [read $in]
close $in
set out [open $name.PRS9a.dat w]
set nBindE 0; set pConE 0; set mConE 0; set aConE 0
set nBindF 0; set pConF 0; set mConF 0; set aConF 0
foreach line [split $inData "\n"] {
  if {$line == ""} break
#  incr frame
  set tm [lindex $line 0]
  incr tm
  set eProt [lindex $line 1]; # if >-1 no eNck protein collision 
  set fProt [lindex $line 2]; # if >-1 no fNck protein collision
  set eMem [lindex $line 3]; # eNck distance to membrane, >-6.0 no collision
  set fMem [lindex $line 4]; # fNck distance to membrane, >-13.0 no collision
  # CD3εδ PRS angles
  set eP182 [lindex $line 5] 
  set eP183 [lindex $line 6] 
  set eP184 [lindex $line 7] 
  # CD3εγ PRS angles
  set fP182 [lindex $line 8]  
  set fP183 [lindex $line 9] 
  set fP184 [lindex $line 10] 
  # No protein contact exlcuding fgr.201
  if {$eProt>-1 && $eProt<24} {incr pConE}
  if {$fProt>-1 && $fProt<24} {incr pConF}
  # No membrane contact
  if {$eMem>-6.0} {incr mConE}
  if {$fMem>-6.0} {incr mConF}
  # 7-19-24 criteria: P182 ψ >50°, P184 ψ >50° permissive for binding.
  # 7-27-24 criteria: P182 ψ >100°, P184 ψ >100° permissive for binding.
  if {$eP182>100 && $eP184>100} {incr aConE}
  if {$fP182>100 && $fP184>100} {incr aConF}
  # Nck can bind if no contacts with other proteins, 
  # no Nck contact with membrane, and no bad angles.
  if {$eProt>-1 && $eProt<24 && $eMem>-6.0 && $eP182>100 && $eP184>100} {incr nBindE}
  if {$fProt>-1 && $fProt<24 && $fMem>-6.0 && $fP182>100 && $fP184>100} {incr nBindF}
  # In 10 ns blocks, each frame is 10 ps.
  if {[expr $tm % 1000]==0} {
    # tcr38.PRS9.dat output in percents 
    # <time (ns)>
    # <eNck no protein contact><fNck no protein contact>
    # <eNck no mem contacts><fNck no mem contacts>
    # <ePRS angles OK><fPRS angles OK>
    # <eNck binding possible><fNck binding possible>
    puts $out [format "%4d %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f" \
        [expr $tm/100] \
        [expr $pConE/10.0] [expr $pConF/10.0] \
        [expr $mConE/10.0] [expr $mConF/10.0] \
        [expr $aConE/10.0] [expr $aConF/10.0] \
        [expr $nBindE/10.0] [expr $nBindF/10.0] ]
    set nBindE 0; set pConE 0; set mConE 0; set aConE 0
    set nBindF 0; set pConF 0; set mConF 0; set aConF 0
  }
}
close $out
}
eof
###########################################################
# Plot Nck binding frequency for Nck binding to chain E and F in tcr38.*.
# Alternate environment variables control which dimer is plotted.
out=tcr38.nckE9.jpg
dimer=CD3εδ
col=0
#
out=tcr38.nckF9.jpg
dimer=CD3εγ 
col=1
#
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out '$out'
# layout: rows, columns
set multiplot layout 3,4 title "{/:Bold Potential Nck binding to $dimer in Tcr38.*}"
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 5
set xtics out nomirror 0,50,200
set ytic out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set key center left textcolor variable samplen -1 font "arial.ttf,12"
# Legend placed in empty panel
set key Left reverse textcolor variable samplen -1 font "arial.ttf,18"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 7 pi -1 ps 1.0
#
noProt=2+$col
noMem=4+$col
favA=6+$col
pot=8+$col
do for [n in "0 2 3 5 7 8 9 10 11"] {
  set title sprintf("tcr38.%s",n)
  fname = sprintf("tcr38.%s.PRS9a.dat",n)
  # Titles in screen coordinates. Last points shifted to expose overlaps.
  xKey = 0.30   
  plot fname using 1:noProt with linespoints ls 1 t "No protein collision" at xKey,0.19,\
  '' using 1:noMem with linespoints ls 2 t "No membrane collision" at xKey, 0.16,\
  '' using 1:favA with linespoints ls 3 t "Favorable PRS angles" at xKey, 0.13,\
  '' using (column(1)+4):pot with linespoints ls 4 t "Potential Nck binding" at xKey, 0.10
  # Offset the last, total binding, to make it visible.
}
quit
eof
###########################################################
# Use R to consolidate Nck binding for all models of tcr38.*.
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# tcr38.*.PRS9.dat has fields above.
# I want the above in tcr38.sum3.dat with the averages over all 9 models. 
files <- Sys.glob("tcr38.*.PRS9a.dat")
nFiles=length(files)
# Matrices of times (10 ns intervals) in rows, binding criteria in columns.
# tcr38.*.PRS9.dat files have 20 or 21 rows. 
Pe<-matrix(,nrow=20,ncol=nFiles)
Pf<-matrix(,nrow=20,ncol=nFiles)
Me<-matrix(,nrow=20,ncol=nFiles)
Mf<-matrix(,nrow=20,ncol=nFiles)
Ae<-matrix(,nrow=20,ncol=nFiles)
Af<-matrix(,nrow=20,ncol=nFiles)
Te<-matrix(,nrow=20,ncol=nFiles)
Tf<-matrix(,nrow=20,ncol=nFiles)
# Fill above matrices with frequencies for each criteria.
for(i in seq(1,nFiles,1)) {
  read.table(files[i])->a
  Pe[,i] <- a[1:20,2]
  Pf[,i] <- a[1:20,3]
  Me[,i] <- a[1:20,4]
  Mf[,i] <- a[1:20,5]
  Ae[,i] <- a[1:20,6]
  Af[,i] <- a[1:20,7]
  Te[,i] <- a[1:20,8]
  Tf[,i] <- a[1:20,9]  
}
# Output has average over all models for binding criteria for each row.
out<-matrix(,nrow=20,ncol=9) 
for(t in seq(1,20,1)) {
  out[t,1]= t*10
  out[t,2] <- round(mean(Pe[t,]),2)
  out[t,3] <- round(mean(Pf[t,]),2)
  out[t,4] <- round(mean(Me[t,]),2)
  out[t,5] <- round(mean(Mf[t,]),2)
  out[t,6] <- round(mean(Ae[t,]),2)
  out[t,7] <- round(mean(Af[t,]),2)
  out[t,8] <- round(mean(Te[t,]),2)
  out[t,9] <- round(mean(Tf[t,]),2)
}
write.table(out,file="tcr38a.sum9.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
# Output: <time><Pe><Pf><Me><Mf><Ae><Af><Te><Tf>
###########################################################
# Plot tcr39b and tcr38b Nck binding averages
# With all criteria for each dimer, four panels.
# Figure S14 A fig
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr3938a.sum9.jpg'
# layout: rows, columns
set multiplot layout 2,2 
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 3.5
set xtics out nomirror 0,50,220 offset 0,0.1
set ytics out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
#set ylabel "{/:Bold Potential Nck binding (%)}" offset 3,0
set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
#
set title "{/:Bold Resting state TCR^Q CD3εδ}" offset 0,-1
plot '/t/tcr39an/tcr39a.sum9.dat' u 1:2 with linespoints ls 1 t "",\
'' using 1:4 with linespoints ls 2 t "" ,\
'' using 1:6 with linespoints ls 3 t "" ,\
'' using 1:8 with linespoints ls 4 t ""  
#
set title "{/:Bold Resting state TCR^Q CD3εγ}" offset 0,-1
plot '/t/tcr39an/tcr39a.sum9.dat' u 1:3 with linespoints ls 1 t "" ,\
'' using 1:5 with linespoints ls 2 t "" ,\
'' using 1:7 with linespoints ls 3 t "" ,\
'' using 1:9 with linespoints ls 4 t ""  
#
set title "{/:Bold Activated TCR-GOF^Q CD3εδ}" offset 0,-1
plot '/t/tcr38an/tcr38a.sum9.dat' u 1:2 with linespoints ls 1 t "" ,\
'' using 1:4 with linespoints ls 2 t "" ,\
'' using 1:6 with linespoints ls 3 t "" ,\
'' using 1:8 with linespoints ls 4 t ""
#
set title "{/:Bold Activated TCR-GOF^Q CD3εγ}" offset 0,-1
set key bottom Left reverse textcolor variable samplen -1 font "arial.ttf,12"
plot '/t/tcr38an/tcr38a.sum9.dat' u 1:3 with linespoints ls 1 t "",\
'' using 1:5 with linespoints ls 2 t "",\
'' using 1:7 with linespoints ls 3 t "",\
'' using 1:9 with linespoints ls 4 t ""  
quit
eof

##
plot '/t/tcr38an/tcr38a.sum9.dat' u 1:3 with linespoints ls 1 t "No protein collision",\
'' using 1:5 with linespoints ls 2 t "No membrane collision",\
'' using 1:7 with linespoints ls 3 t "Favorable PRS angles",\
'' using 1:9 with linespoints ls 4 t "Potential Nck binding" 
