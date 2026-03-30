# anEdt20.tcl: Test Lck binding to edt 200 ns MD (Fig 7A,C).
# Output edtLck1.$mod.dat
# Derived from anEdt7.tcl 
# vmd -dispdev none on Eire.
# Contacts are heavy atoms less than or equal to 1.0 angstroms apart. 
set rmdsMax 3.69;  # Limit for Lck substrate fit.
# ^ Not used here
cd /t/tcr16/edtAn
###########################################################
# Lck substrate atomselection
set lckpep [mol new /home/jcannon/lck/Lck-Peptide.pdb]
set lck [atomselect $lckpep "chain A"]
set lcksub [atomselect $lckpep all]
set lckH [atomselect $lckpep "chain A and mass>2"]; # Lck heavy atoms
# Choose residues for fitting
set t 4; # The 9 residues around pTyr for fitting.
set w 7; # How far away from pTyr are Lck collisions allowed? 
set c [expr 285- $t]; # first Lck substrate residue to fit
set d [expr 285+ $t]; # last Lck substrate residue to fit
set pSel [atomselect $lckpep "name CA and chain B and resid $c to $d"]
###########################################################
# ITAM list with {chain resid}
set itamsAll {{A 72} {A 83} {A 111} {A 123} {A 142} {A 153} 
{B 72} {B 83} {B 111} {B 123} {B 142} {B 153}  
{D 149} {D 160} {E 188} {E 199} {F 188} {F 199} {G 160} {G 171} }
set itams { {D 149} {D 160} {E 188} {E 199} }
###########################################################
# Load edt 5.4-205 ns MD trajectories for 11 edt models.
#set edtModels "84 171 30 115 35 90 13 41 32 153 198"
set edtModels "35"
foreach mod $edtModels {
set out [open edtLck1.$mod.dat w]
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
set lastframe [molinfo $tcr get numframes]
# CA of first chain E CT residue Ala157 used to judge membrane location. 
set Mem [atomselect $tcr "name CA and chain E and resid 157"]
###########################################################
# ITAM atomselection's
set i 0
foreach s $itams {
set chn [lindex $s 0]; # ITAM chain
set n [lindex $s 1];   # ITAM resid
set a [expr $n- $t]; # first ITAM residue to fit
set b [expr $n+ $t]; # last ITAM residue to fit
set iSel($i) [atomselect $tcr "chain $chn and name CA and resid $a to $b"]
# Check proximal iNot (same as ITAM chain) or other protein colissions jNot.
set e [expr $n- $w];set g [expr $n+ $w]
set iNot($i) [atomselect $tcr "chain $chn and mass>2 and not (resid $e to $g)"]
set jNot($i) [atomselect $tcr "protein and mass>2 and not chain $chn"]
incr i
}
# 
for {set f 0} {$f<$lastframe} {incr f} {
  # Start output line with frame#
  puts -nonewline $out [format "%5d" $f]
  $Mem frame $f
  set memZ [$Mem get z]
  for {set k 0} {$k<$i} {incr k} {
    $iSel($k) frame $f; $iNot($k) frame $f; $jNot($k) frame $f
    # Fitting: fit substrate to ITAM and move lcksub
    set M [measure fit $pSel $iSel($k)]
    $lcksub move $M
    set rms [measure rmsd $pSel $iSel($k)]
    # Check proximal, other proteins, and membrane Lck collisions
    set conSub1 [measure contacts 1.0 $lckH $iNot($k)]
    set prox [llen [lindex $conSub1 1]]
    set conSub2 [measure contacts 1.0 $lckH $jNot($k)]
    set other [llen [lindex $conSub2 1]]
    set conSub3 [measure contacts 1.0 $lckH $Mem]
    set memC [llen [lindex $conSub3 1]]
    # ^ That memC above does not count membrane contacts. 
    # Instead, use dMem below. It must be >-6.0 to not contact membrane.
    # Get Lck to membrane distance
    set loc [measure minmax $lcksub]
    set minZ [lindex $loc 0 2]
    set dMem [expr $minZ- $memZ]
    # Output for each ITAM Tyr: <rmsFit><Lck distance to membrane><proximal contacts>
    # <$other protein contacts><membrane contacts>
    puts -nonewline $out [format "%5.2f %5.2f %4d%4d%4d | " $rms $dMem $prox $other $memC ]
  }
  puts $out " "; # To finish output line
  if {[expr $f % 10]==0} {puts "$f edt.$mod"}; # progress report
}  
close $out
mol delete $tcr
for {set k 0} {$k<$i} {incr k} {
    $iSel($k) delete; $iNot($k) delete; $jNot($k) delete
}
}
###########################################################
# Analysis for frequency of Lck binding to edt frames.
# Report individual criteria as well as cumulative criteria for each of four ITAM Tyr.
# Set d to D149, D160, E188, or E199 for each of four ITAM Tyr. 
tclsh<<"eof"
set models {84 171 30 115 35 90 13 41 32 153 198} 
foreach mod $models {
puts $mod
set in [open edtLck1.$mod.dat r]
set inData [read $in]
close $in
# Set d to ITAM Tyr chain, resid
foreach d "D149 D160 E188 E199" {
set out [open edt$d.$mod.dat w]
# Example output of edtLck1.$mod.dat. 
# <frame><rmsFit><Lck distance to membrane><proximal contacts>
    # <other protein contacts><membrane contacts>
# 0 2.61  10.40  3  4  0 |  3.13  31.94  0  0  0 |  2.84  -8.23 14  0  0 |  2.70  -3.78  5  0  0 |
# 1 2.84   4.99 10  7  0 |  3.25  28.57  0  0  0 |  2.95  -8.24 17  0  0 |  2.79  -4.03  5  0  0 |
# Data for each ITAM Tyr in this order (chain resid): D 149; D 160; E 188; E 199.
set pconT 0; set oconT 0; set memT 0; set rmsT 0; 
set bin4T 0; set bin3T 0; set bin2T 0
foreach line [split $inData "\n"] {
  if {$line == ""} break
  set a [llength $line]
  # Check for the correct number of fields on each line.
  if {$a != 25} {puts "ERROR: too few fields";puts [lindex $line 0]; exit}
  set tm [lindex $line 0]
  incr tm
  if {[string equal $d D149]} {
    set rms [lindex $line 1]
    set dis [lindex $line 2]
    set pcon [lindex $line 3]
    set ocon [lindex $line 4]
    set mem [lindex $line 5]
  }
  if {[string equal $d D160]} {
    # chain D Tyr160
    set rms [lindex $line 7]
    set dis [lindex $line 8]
    set pcon [lindex $line 9]
    set ocon [lindex $line 10]
    set mem [lindex $line 11]
  }
  if {[string equal $d E188]} {
    # chain e Tyr188
    set rms [lindex $line 13]
    set dis [lindex $line 14]
    set pcon [lindex $line 15]
    set ocon [lindex $line 16]
    set mem [lindex $line 17]
  }
  if {[string equal $d E199]} {
    # chain e Tyr199
    set rms [lindex $line 19]
    set dis [lindex $line 20]
    set pcon [lindex $line 21]
    set ocon [lindex $line 22]
    set mem [lindex $line 23]
  }
  if {$rms<3.69} {incr rmsT}
  if {$pcon==0} {incr pconT}
  if {$ocon==0} {incr oconT}
  # membrane contacts are not used, use Lck distance to membrane.
#  if {$mem==0} {incr memT}
  if {$dis>-6} {incr memT}
  # 4 criteria: Lck can bind if rms<3.69 and no protein or membrane contacts.
  if {$rms<3.69 && $pcon==0 && $ocon==0 && $dis>-6} {incr bin4T}
  # 3 criteia: Lck can bind if rms1<3.69 and no other protein or membrane contacts.
  if {$rms<3.69 && $ocon==0 && $dis>-6} {incr bin3T}
  # 2 criteria: Lck can bind if rms1<3.69 and no proximal contacts.
  if {$rms<3.69 && $pcon==0} {incr bin2T} 
  # In 10 ns blocks, each frame is 10 ps.
  if {[expr $tm % 1000]==0} {
    # Output in percents 
    # <time (ns)><no proximal contact><no other p contact><no mem contacts>\
    # <RMSD OK><binding possible with proximal><binding possible without proximal>
    # <RMSD and proximal OK>
    puts $out [format "%4d %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f" \
        [expr $tm/100] [expr $pconT/10.0] [expr $oconT/10.0] \
        [expr $memT/10.0] [expr $rmsT/10.0] [expr $bin4T/10.0] \
        [expr $bin3T/10.0] [expr $bin2T/10.0]]
    set pconT 0; set oconT 0; set memT 0; set rmsT 0; 
    set bin4T 0; set bin3T 0; set bin2T 0
  }
}
close $out
}
}
eof
###########################################################
# Plot Lck binding frequency for eleven edt models in multiplot
# edt models in cluster prevalence order:
models="84 171 30 115 35 90 13 41 32 153 198"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
# Set c to 3, 4 for # of Lck binding criteria
c=4
# Set d to D149, D160, E188, or E199 for Tyr chain, resid.
d="E199"
outName = sprintf("edt%d.%s.jpg",c,d)
set out outName
# layout: rows, columns
set multiplot layout 3,4 title sprintf("edt%d.%s",c,d)
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 5
set xtics out nomirror 0,50,220
set ytic out nomirror offset 1,0
set xrange [0:210]
set yrange [0:100]
# Legend placed in lower right empty panel.
set key Left reverse textcolor variable samplen -1 font "arial.ttf,18"
# These have a variety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
#
do for [n in "$models"] {
  set title sprintf("edt.%s",n) offset 0,-1
  fname = sprintf("edt%s.%s.dat",d,n)
  # Titles in screen coordinates. Last points shifted to expose overlaps.
  # Offset the last, total binding to make it visible.
  if (c==4) {
  # 4 criteria
  plot fname using 1:2 with linespoints ls 1 t "No proximal collision" at 0.77,0.22,\
  '' using 1:3 with linespoints ls 5 t "No other collision" at 0.77, 0.19,\
  '' using 1:4 with linespoints ls 2 t "No membrane collision" at 0.77, 0.16,\
  '' using 1:5 with linespoints ls 3 t "Favorable RMSD" at 0.77, 0.13,\
  '' using (column(1)+4):6 with linespoints ls 4 t "Potential Lck binding" at 0.77, 0.10  
  }
  if (c==3) {
  # 3 criteria   
  plot fname using 1:3 with linespoints ls 5 t "No other collision" at 0.77, 0.19,\
  '' using 1:4 with linespoints ls 2 t "No membrane collision" at 0.77, 0.16,\
  '' using 1:5 with linespoints ls 3 t "Favorable RMSD" at 0.77, 0.13,\
  '' using (column(1)+4):7 with linespoints ls 4 t "Potential Lck binding" at 0.77, 0.10
  }
}
quit
eof
###########################################################
# Use R to consolidate Lck binding for all models of edt.
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# Files edt$d.$mod.dat (d=D149, D160, E188, E199; mod=model#) have
# <time (ns)><no proximal contact><no other p contact><no mem contacts>\
# <RMSD OK><binding possible with proximal><binding possible without proximal>
# Output a separate summary for each ITAM Tyr.
# Set d to D149, D160, E188, or E199 for Tyr chain, resid.
itams=c("D149", "D160", "E188", "E199")
for(d in itams) {
# I want the above in edtSum.dat with the averages over all 11 models. 
models=c("84","171","30","115","35","90","13","41","32","153","198")
nMod=length(models)
nTime=20
# Matrices of times in rows, models in columns
P<-matrix(,nrow=nTime,ncol=nMod)
O<-matrix(,nrow=nTime,ncol=nMod)
M<-matrix(,nrow=nTime,ncol=nMod)
R<-matrix(,nrow=nTime,ncol=nMod)
T<-matrix(,nrow=nTime,ncol=nMod)
Q<-matrix(,nrow=nTime,ncol=nMod)
F<-matrix(,nrow=nTime,ncol=nMod)
# Fill above matrices with frequencies for each criteria.
for(i in seq(1,nMod,1)) {
  read.table(sprintf("edt%s.%s.dat",d,models[i]))->a
  P[,i] <- a[,2]
  O[,i] <- a[,3]
  M[,i] <- a[,4]
  R[,i] <- a[,5]
  T[,i] <- a[,6]
  Q[,i] <- a[,7]
  F[,i] <- a[,8]
}
# Output has average over all models for binding criteria for each row.
# Each column has a different criterion.
out<-matrix(,nrow=nTime,ncol=8) 
for(t in seq(1,nTime,1)) {
  out[t,1]= t*10
  out[t,2] <- round(mean(P[t,]),2)
  out[t,3] <- round(mean(O[t,]),2)
  out[t,4] <- round(mean(M[t,]),2)
  out[t,5] <- round(mean(R[t,]),2)
  out[t,6] <- round(mean(T[t,]),2)
  out[t,7] <- round(mean(Q[t,]),2)
  out[t,8] <- round(mean(F[t,]),2)
}
write.table(out,file=sprintf("edtLck%s.dat",d),quote=FALSE,row.names=FALSE,col.names=FALSE)
} # foreach itams
###########################################################
# Plot Lck binding averages for edt ensemble using four criteria
# Fig 7A
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edt.lck1.jpg'
# layout: rows, columns
set multiplot layout 2,2 title "{/:Bold Potential Lck binding to CD3εδ}"
set border 3
set lmargin 4
set rmargin 2
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
#set ylabel "{/:Bold Potential Lck binding (%)}" offset 3,0
set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
# These have a variety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
xKey=0.8; # x coordinate for legend
do for [d in "D149 D160 E188 E199"] {
  if (d eq "D149") {set label sprintf("{/:Bold CD3δ Tyr149}",d) at graph 0.05,1.0}
  if (d eq "D160") {set label sprintf("{/:Bold CD3δ Tyr160}",d) at graph 0.05,1.0}
  if (d eq "E188") {set label sprintf("{/:Bold CD3ε Tyr188}",d) at graph 0.05,1.0}
  if (d eq "E199") {set label sprintf("{/:Bold CD3ε Tyr199}",d) at graph 0.05,1.0} 
  fname = sprintf("edtLck%s.dat",d)
  plot fname using 1:2 with linespoints ls 1 t "No proximal collision" at xKey,0.48,\
  '' using 1:3 with linespoints ls 5 t "No other collision" at xKey, 0.46,\
  '' using 1:4 with linespoints ls 2 t "No membrane collision" at xKey, 0.44,\
  '' using 1:5 with linespoints ls 3 t "Favorable RMSD" at xKey, 0.42,\
  '' using (column(1)+4):6 with linespoints ls 4 t "Potential Lck binding" at xKey, 0.40 
  unset label 
}
quit
eof
###########################################################
# Plot Lck binding averages for edt ensemble using two criteria
# Fig 7C
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edt.lck6.jpg'
# layout: rows, columns
set multiplot layout 2,2 title "{/:Bold Potential Lck binding to free CD3εδ}"
set border 3
set lmargin 4
set rmargin 2
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
#set ylabel "{/:Bold Potential Lck binding (%)}" offset 3,0
set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
# These have a variety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
xKey=0.8; # x coordinate for legend
do for [d in "D149 D160 E188 E199"] {
  if (d eq "D149") {set label sprintf("{/:Bold CD3δ Tyr149}",d) at graph 0.05,1.0}
  if (d eq "D160") {set label sprintf("{/:Bold CD3δ Tyr160}",d) at graph 0.05,1.0}
  if (d eq "E188") {set label sprintf("{/:Bold CD3ε Tyr188}",d) at graph 0.05,1.0}
  if (d eq "E199") {set label sprintf("{/:Bold CD3ε Tyr199}",d) at graph 0.05,1.0} 
  fname = sprintf("edtLck%s.dat",d)
  plot fname using 1:2 with linespoints ls 1 t "No proximal collision" at xKey,0.44,\
  '' using 1:5 with linespoints ls 3 t "Favorable RMSD" at xKey, 0.42,\
  '' using (column(1)+4):8 with linespoints ls 4 t "Potential Lck binding" at xKey, 0.40 
  unset label 
}
quit
eof
# Use R to get statistics on above data
read.table("edtLckD149.dat")->a
mean(a[,5])
mean(a[,2])
mean(a[,8])
read.table("edtLckD160.dat")->a
mean(a[,5])
mean(a[,2])
mean(a[,8])
read.table("edtLckE188.dat")->a
mean(a[,5])
mean(a[,2])
mean(a[,8])
read.table("edtLckE199.dat")->a
mean(a[,5])
mean(a[,2])
mean(a[,8])
