# anFgt14.tcl: Nck binding to fgt (S6 fig panel B).
# Use vmd -dispdev none on Eire.
# Check CD3 PRS dihedral angles for Nck binding criteria
# and combine with collision data to calculate Nck binding.
# Use the 7-27-24 criteria: P182 ψ <100°, P184 ψ <100°
# This adds to the anfgt7.tcl output, fgtPRS10.*.dat 
# This is derived from anFgt14.tcl
# Output fgtPRS16.$mod.dat 
# vmd -dispdev none on Eire.
cd /t/tcr16/fgtAn
# Using wrapped angles, these limits prevent Nck binding
# P182 ψ <50°, V183 φ<240°, P184 ψ <50° 
proc wrapPhi {a} {
  return [expr $a<0?$a+360:$a]
}
proc wrapPsi {a} {
  return [expr $a<-120?$a+360:$a]
}
###########################################################
# Load fgt 5.4-205 ns MD trajectories for 14 fgt models.
set fgtModels "185 121 191 180 137 75 14 153 93 28 109 122 42 45"
foreach mod $fgtModels {
set out [open fgtPRS16.$mod.dat w]
# Use stripped topology with just 2448 protein atoms that were saved in trajectory.
# Made by anFgt1.tcl
set tcr [mol new /t/tcr16/fgt/fgt.prot.psf]
# The NPT MD started with fgt.*.eq2.
for {set i 2} {$i<=55} {incr i} {
  set mdFile /t/tcr16/fgt/fgt.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
# Use fgt resid from fgt.prot.psf, which is different from fgt.*.top
set P182 [atomselect $tcr "chain F and name CA and resid 182"]
set P184 [atomselect $tcr "chain F and name CA and resid 184"]
set lastframe [molinfo $tcr get numframes]
for {set f 0} {$f<$lastframe} {incr f} {
  set badA 0; set bad182 0; set bad184 0
  $P182 frame $f; $P184 frame $f
  # Use the 7-27-24 criteria: P182 ψ <100°, P184 ψ <100°
  set P182psi [$P182 get psi] 
  if {[wrapPsi $P182psi] < 100} {incr badA; incr bad182}
  set P184psi [$P184 get psi]
  if {[wrapPsi $P184psi] < 100} {incr badA; incr bad184}
  # Output: <frame><total bad angles><bad 182><bad184>
  puts $out [format "%6d%4d%4d%4d" $f $badA $bad182 $bad184]
}
close $out
mol delete $tcr
}
###########################################################
# Combine collision in fgtPRS10.*.dat with bad angles in fgtPRS14.*.dat
for mod in 185 121 191 180 137 75 14 153 93 28 109 122 42 45; do
paste fgtPRS10.$mod.dat fgtPRS16.$mod.dat >temp
# Save <frame><nRef><RMSD fit><membrane dist><bad angles><bad 182><bad 184>
awk '{printf("%6d%4d%7.2f%7.2f%4d%4d%4d\n",$1,$2,$3,$4,$6,$7,$8)}' \
temp>fgtPRS16.$mod.dat 
done
rm -f temp
# Now fgtPRS16.$mod.dat has
# <frame><nRef><RMSD fit><eNck membrane><bad angles>
###########################################################
# Analysis for frequency of Nck binding in 10 ns blocks. 
tclsh<<"eof"
set models {185 121 191 180 137 75 14 153 93 28 109 122 42 45} 
foreach mod $models {
puts $mod
set in [open fgtPRS16.$mod.dat r]
set inData [read $in]
close $in
set out [open fgtPRS17.$mod.dat w]
# Example output of fgtPRS16.$mod.dat. 
#      0   0   2.33  21.55   1   0   1   0
set nBind 0; set pCon 0; set mCon 0; set aCon 0
foreach line [split $inData "\n"] {
  if {$line == ""} break
  incr frame
  set tm [lindex $line 0]
  incr tm
  # Get Nck-PRS ref index
  set nRef [lindex $line 1]; #  Nck-PRS ref index 
  incr r($nRef); # track how many times each Nck-PRS ref used.
  set RMSD [lindex $line 2]; # Disregard the PRS RMSD
  set dMem [lindex $line 3]; # Distance to membrane
  set badA [lindex $line 4]; # number of bad PRS angles
  set bad182 [lindex $line 5]; # bad 182
  set bad184 [lindex $line 6]; # bad 184
  # No protein contact exlcuding fgr.201
  if {$nRef>-1 && $nRef<24} {incr pCon}
  # No membrane contact
  if {$dMem>-13.0} {incr mCon}
  # Just look at P182 and P184
  if {$bad182==0 && $bad184==0} {incr aCon}
  # Nck can bind if no contacts with other proteins, 
  # no Nck contact with membrane, and no bad angles.
  if {$nRef>-1 && $nRef<24 && $dMem>-6.0 && $bad182==0 && $bad184==0} {incr nBind}
  # In 10 ns blocks, each frame is 10 ps.
  if {[expr $tm % 1000]==0} {
    # Output in percents 
    # <time (ns)><no protein contact><no mem contacts><angles OK><binding possible>
    puts $out [format "%4d %4.3f %4.3f %4.3f %4.3f" \
        [expr $tm/100] [expr $pCon/10.0] [expr $mCon/10.0] \
        [expr $aCon/10.0] [expr $nBind/10.0]]
    set nBind 0; set pCon 0; set mCon 0; set aCon 0
  }
}
close $out
}
# Output how many times each Nck reference bound.
# This was useful in the ordering of references in anfgt7.tcl.
parray r
puts "Total frames=$frame"
eof
###########################################################
# Plot Nck binding frequency for eleven fgt models in multiplot
# Figure S5B
# fgt models in cluster prevalence order:
models="185 121 191 180 137 75 14 153 93 28 109 122 42 45"
gnuplot<<eof
#set term jpeg font "arial.ttf,18" size 1280,960
set term jpeg font "arial.ttf,18" size 1000,750
set out 'fgt.nck16.jpg'
# layout: rows, columns
#set multiplot layout 4,4 title "{/:Bold Potential Nck binding to CD3εγ}"
set multiplot layout 4,4
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 3.5
set tmargin 1
set xtics out nomirror 0,100,220
set xtics add ("" 50,"" 150)
set ytic out nomirror offset 1,0
set xrange [0:210]
set yrange [0:100]
set grid y
# Legend placed in lower right 
set key Left reverse textcolor variable samplen -1 font "arial.ttf,16"
# These have a vairety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
#
L=0.6
do for [n in "$models"] {
  set title sprintf("fgt.%s",n) offset 0,-1
  fname = sprintf("fgtPRS17.%s.dat",n)
  # Columns in percent <time (ns)><no protein contact><no contacts><Nck bound>
  # Titles in screen coordinates 
  plot fname using 1:2 with linespoints ls 1 t "No protein collision" at L,0.19,\
  '' using 1:3 with linespoints ls 2 t "No membrane collision" at L, 0.16,\
  '' using 1:4 with linespoints ls 3 t "Favorable PRS angles" at L, 0.13,\
  '' using (column(1)+4):5 with linespoints ls 4 t "Potential Nck binding" at L, 0.10
}
quit
eof
###########################################################
# Use R to consolidate Nck binding for all models of fgt.
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# fgtPRS17.*.dat has
# <time (ns)><no protein contact><no mem contacts><angles OK><binding possible>
# I want the above in fgtSum.dat with the averages over all 14 models. 
models=c("185","121","191","180","137","75","14","153","93","28","109","122","42","45")
# Matrices of times in rows, models in columns
# The fgt trajectories wnet to 210 ns, so there are 21 rows.
P<-matrix(,nrow=21,ncol=14)
M<-matrix(,nrow=21,ncol=14)
A<-matrix(,nrow=21,ncol=14)
T<-matrix(,nrow=21,ncol=14)
# Fill above matrices with frequencies for each criteria.
for(i in seq(1,14,1)) {
  read.table(sprintf("fgtPRS17.%s.dat",models[i]))->a
  # The fgt trajectories wnet to 210 ns, so there are 21 rows.
  P[,i] <- a[,2]
  M[,i] <- a[,3]
  A[,i] <- a[,4]
  T[,i] <- a[,5]
}
# Output has average over all models for binding criteria for each row.
out<-matrix(,nrow=21,ncol=5) 
for(t in seq(1,21,1)) {
  out[t,1]= t*10
  out[t,2] <- round(mean(P[t,]),2)
  out[t,3] <- round(mean(M[t,]),2)
  out[t,4] <- round(mean(A[t,]),2)
  out[t,5] <- round(mean(T[t,]),2)
}
write.table(out,file="fgtSum.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
###########################################################
# Plot fgt Nck binding averages,
gnuplot<<eof
set term jpeg font "arial.ttf,18" 
set out 'fgt.nck16b.jpg'
set title "{/:Bold Potential Nck binding to CD3εγ in fgt}"
set border 3
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
set ylabel "{/:Bold Potential Nck binding (%)}" offset 3,0
set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
# These have a vairety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
plot '/t/tcr16/fgtAn/fgtSum.dat' u 1:2 with linespoints ls 1 t "No protein collision" at 0.25,0.29,\
'' using 1:3 with linespoints ls 2 t "No membrane collision" at 0.25,0.26,\
'' using 1:4 with linespoints ls 3 t "Favorable PRS angles" at 0.25,0.23,\
'' using (column(1)+4):5 with linespoints ls 4 t "Potential Nck binding" at 0.25,0.20 
quit
eof
