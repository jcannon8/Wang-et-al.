# anEdt16.tcl: Check CD3 PRS dihedral angles for Nck binding criteria
# and combine with collision data to calculate Nck binding.
# Use the 7-27-24 criteria: P182 ψ <100°, P184 ψ <100°
# This is derived from anEdt14.tcl
# This uses anEdt7.tcl collision output, edtPRS10.*.dat 
# This outputs edtPRS16.$mod.dat 
# vmd -dispdev none on Eire.
cd /t/tcr16/edtAn
# Using wrapped angles,  
proc wrapPhi {a} {
  return [expr $a<0?$a+360:$a]
}
proc wrapPsi {a} {
  return [expr $a<-120?$a+360:$a]
}
###########################################################
# Load edt 5.4-205 ns MD trajectories for 11 edt models.
set edtModels "84 171 30 115 35 90 13 41 32 153 198"
foreach mod $edtModels {
set out [open edtPRS16.$mod.dat w]
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
# Use edt resid from edt.prot.psf, which is different from edt.*.top
set P182 [atomselect $tcr "chain E and name CA and resid 182"]
set P184 [atomselect $tcr "chain E and name CA and resid 184"]
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
# Combine collision in edtPRS10.*.dat with bad angles in edtPRS16.*.dat
for mod in 84 171 30 115 35 90 13 41 32 153 198; do
paste edtPRS10.$mod.dat edtPRS16.$mod.dat >temp
# Save <frame><nRef><RMSD fit><membrane dist><bad angles><bad 182><bad 184>
awk '{printf("%6d%4d%7.2f%7.2f%4d%4d%4d\n",$2,$3,$4,$5,$7,$8,$9)}' \
temp>edtPRS16.$mod.dat 
done
rm -f temp
# Now edtPRS16.$mod.dat has
# <frame><nRef><RMSD fit><membrane dist><bad angles><bad 182><bad184>
cat edtPRS16.*.dat > edtPRS16.sum.dat
###########################################################
# Analysis for frequency of Nck binding in 10 ns blocks.
# This reports individual criteria as well as cumulative criteria. 
tclsh<<"eof"
set models {84 171 30 115 35 90 13 41 32 153 198} 
foreach mod $models {
puts $mod
set in [open edtPRS16.$mod.dat r]
set inData [read $in]
close $in
set out [open edtPRS17.$mod.dat w]
# Example output of edtPRS16.$mod.dat. 
#      0   0   2.33  21.55   1   0   1   0
# <frame><nRef><RMSD fit><membrane dist><bad angles><bad 182><bad184>
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
  if {$dMem>-6.0} {incr mCon}
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
# This was useful in the ordering of references in anEdt7.tcl.
parray r
puts "Total frames=$frame"
eof
###########################################################
# Plot Nck binding frequency for eleven edt models in multiplot
# Figure S5A
# edt models in cluster prevalence order:
models="84 171 30 115 35 90 13 41 32 153 198"
gnuplot<<eof
#set term jpeg font "arial.ttf,18" size 1280,960
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edt.nck16.jpg'
# layout: rows, columns
#set multiplot layout 3,4 title "{/:Bold Potential Nck binding to CD3εδ}"
set multiplot layout 3,4
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 3.5
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
L=0.75
do for [n in "$models"] {
  set title sprintf("edt.%s",n) offset 0,-1
  fname = sprintf("edtPRS17.%s.dat",n)
  # Titles in screen coordinates. Last points shifted to expose overlaps.   
  plot fname using 1:2 with linespoints ls 1 t "No protein collision" at L,0.19,\
  '' using 1:3 with linespoints ls 2 t "No membrane collision" at L, 0.16,\
  '' using 1:4 with linespoints ls 3 t "Favorable PRS angles" at L, 0.13,\
  '' using (column(1)+4):5 with linespoints ls 4 t "Potential Nck binding" at L, 0.10
  # Offset the last, total binding to make it visible.
}
quit
eof
###########################################################
# Plot RMS vs badA
models="84 171 30 115 35 90 13 41 32 153 198"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out 'edt.nck18.jpg'
# layout: rows, columns
set multiplot layout 3,4 title "{/:Bold RMS vs badA in CD3εδ}"
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 5
set xtics out nomirror 0,1,3
set ytic out nomirror offset 1,0
set grid xtics ytics
set xrange [-0.5:3.5]
set yrange [0:5]
# Legend placed in lower right 
set key Left reverse textcolor variable samplen -1 font "arial.ttf,18"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "dark-cyan" lt 1 lw 2 pt 7 pi -1 ps 1.0
#
do for [n in "$models"] {
  set title sprintf("edt.%s",n) offset 0,-1
  fname = sprintf("edtPRS16.%s.dat",n) 
  # Plot average RMS y-value for each badA x-value
  plot fname using 5:3 with linespoints s unique ls 2 t "" 
}
quit
eof
###########################################################
# Use R to consolidate Nck binding for all models of edt.
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# edtPRS17.*.dat has
# <time (ns)><no protein contact><no mem contacts><angles OK><binding possible>
# I want the above in edtSum.dat with the averages over all 11 models. 
models=c("84","171","30","115","35","90","13","41","32","153","198")
# Matrices of times in rows, models in columns
P<-matrix(,nrow=20,ncol=11)
M<-matrix(,nrow=20,ncol=11)
A<-matrix(,nrow=20,ncol=11)
T<-matrix(,nrow=20,ncol=11)
# Fill above matrices with frequencies for each criteria.
for(i in seq(1,11,1)) {
  read.table(sprintf("edtPRS17.%s.dat",models[i]))->a
  P[,i] <- a[,2]
  M[,i] <- a[,3]
  A[,i] <- a[,4]
  T[,i] <- a[,5]
}
# Output has average over all models for binding criteria for each row.
out<-matrix(,nrow=20,ncol=5) 
for(t in seq(1,20,1)) {
  out[t,1]= t*10
  out[t,2] <- round(mean(P[t,]),2)
  out[t,3] <- round(mean(M[t,]),2)
  out[t,4] <- round(mean(A[t,]),2)
  out[t,5] <- round(mean(T[t,]),2)
}
write.table(out,file="edtSum.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
###########################################################
cd /t/tcr16/edtAn
# Plot edt Nck binding averages for ensemble
gnuplot<<eof
#set term jpeg font "arial.ttf,18" size 1000,750
set term jpeg font "arial.ttf,18" size 500,375
set out 'edt.nck16b.jpg'
set title "{/:Bold CD3εδ}"
set border 3
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
#set ylabel "{/:Bold Potential Nck binding (%)}" offset 3,0
set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
# These have a variety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
#
L=0.25
plot '/t/tcr16/edtAn/edtSum.dat' u 1:2 with linespoints ls 1 t "" ,\
'' using 1:3 with linespoints ls 2 t "" ,\
'' using 1:4 with linespoints ls 3 t "" ,\
'' using (column(1)+4):5 with linespoints ls 4 t "" 
quit
eof
#
plot '/t/tcr16/edtAn/edtSum.dat' u 1:2 with linespoints ls 1 t "No protein collision" at L,0.29,\
'' using 1:3 with linespoints ls 2 t "No membrane collision" at L,0.26,\
'' using 1:4 with linespoints ls 3 t "Favorable PRS angles" at L,0.23,\
'' using (column(1)+4):5 with linespoints ls 4 t "Potential Nck binding" at L,0.20 
###########################################################
# Fig 5 four panel plot
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edt.nck16c.jpg'
# layout: rows, columns
set multiplot layout 2,2 
set border 3
set rmargin 1
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set xrange [0:220]
set yrange [0:100]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
#set ylabel "{/:Bold Potential Nck binding (%)}" offset 3,0
set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
# These have a variety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
#
set title "{/:Bold CD3εδ}" offset 0,-1
plot '/t/tcr16/edtAn/edtSum.dat' u 1:2 with linespoints ls 1 t "" ,\
'' using 1:3 with linespoints ls 2 t "" ,\
'' using 1:4 with linespoints ls 3 t "" ,\
'' using (column(1)+4):5 with linespoints ls 4 t "" 
#
set title "{/:Bold CD3εγ}" offset 0,-1
plot '/t/tcr16/fgtAn/fgtSum.dat' u 1:2 with linespoints ls 1 t "",\
'' using 1:3 with linespoints ls 2 t "" ,\
'' using 1:4 with linespoints ls 3 t "" ,\
'' using (column(1)+4):5 with linespoints ls 4 t "" 
#
set xrange [0:110]
set title "{/:Bold CD3εδ^S}" offset 0,-1
L=0.25
plot '/t/tcr16/eds/edsSum.dat' u 1:2 with linespoints ls 1 t "No protein collision" at L,0.22,\
'' using 1:3 with linespoints ls 2 t "No membrane collision" at L,0.19,\
'' using 1:4 with linespoints ls 3 t "Favorable PRS angles" at L,0.16,\
'' using 1:5 with linespoints ls 4 t "Potential Nck binding" at L,0.13 
#
set title "{/:Bold CD3εγ^S}" offset 0,-1
plot '/t/tcr16/fgs/fgsSum.dat' u 1:2 with linespoints ls 1 t "" ,\
'' using 1:3 with linespoints ls 2 t "",\
'' using 1:4 with linespoints ls 3 t "",\
'' using 1:5 with linespoints ls 4 t "" 
quit
eof
