# tcr38lck3.tcl: Plot Lck binding to tcr38b, tcr39b (Fig 16A,B).
# Using data from tcr38lck2.tcl
# Example of tcr38lck2.dat output:
#   tcr38.0.00   8 A  153  3.12    9.57    9.45
# <model><MD frame eq><chain><resid><RMSD><to mem mean><to mem>
cd /t/tcr38an
tclsh<<"eof"
set in [open tcr38lck2.dat r]
set outA [open tcr38lck3a.dat w]; # Output for Lck binding counts
set outB [open tcr38lck3b.dat w]; # Output for Lck to membrane distances
set outC [open tcr38lck3c.dat w]; # Output for ITAM ID's
set inData [read $in]
close $in
foreach line [split $inData "\n"] {
  if {$line == ""} break
  set mod [lindex $line 0]
  set f [lindex $line 1]; # frame, eq#
  set chn [lindex $line 2]; # ITAM chain
  set res [lindex $line 3]; # ITAM resid
  set rms [lindex $line 4]; # RMSD of Lck-peptide fit
  set meanMem [lindex $line 5]; # Lck to mean membrane
  set LckMem [lindex $line 6]; # Lck to membrane
  #
  incr b($f);  # Increment binding in this frame
  lappend m($f) $LckMem;  # Add Lck to membrane distance to list.
  set itam [format "%s %s" $chn $res];  # construct ITAM ID
  if {$f<=45} {incr idA($itam)};  # Increment binding to this ITAM eq4-45
  if {$f>45} {incr idB($itam)};  # Increment binding to this ITAM eq46-85 
}
# Output collected data
for {set i 4} {$i<=85} {incr i} {
  # 62 models per frame, 20 ITAM Tyr per model.
  puts $outA [format "%3d %3d %4.3f" \
    [expr $i*5] $b($i) [expr 100.0*$b($i)/(20*62)]]; # <time><# Lck bindings><% Lck bound>
  puts $outB $m($i); # List of Lck to membrane distances for this frame
}
# Per ITAM bindings
set itams {{A 72} {A 83} {A 111} {A 123} {A 142} {A 153} 
{B 72} {B 83} {B 111} {B 123} {B 142} {B 153}  
{D 149} {D 160} {E 188} {E 199} {F 188} {F 199} {G 160} {G 171} }
foreach s $itams {
  set chn [lindex $s 0]
  set res [lindex $s 1]
  set itam [format "%s %s" $chn $res];  # construct ITAM ID
  # 62 models, 40 frames
  puts $outC [format "%7s %4.3f %4.3f" \
    $itam [expr 100.0*$idA($itam)/(62*40)] [expr 100.0*$idB($itam)/(62*40)]]
}
close $outA; close $outB; close $outC
# maximum Lck bindings per frame = 39
eof
###########################################################
# Plot tcr39b and tcr38b Lck binding averages for Figure TcrLck2 A
# Fig 16A
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 500,375
set out 'tcr.lck3.jpg'
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 4
set xtics out nomirror 0,100,430 add ("" 50,"" 150,"" 250,"" 350) offset 0,0.2
set ytics out nomirror 0,1,4 add ("" 0.5,"" 1.5,"" 2.5,"" 3.5) offset 1,0
set xrange [0:430]
set yrange [0:4]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
set ylabel "{/:Bold Potential Lck binding (%)}" offset 2,-0.5
set key Left top reverse textcolor variable samplen -1 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
#
plot '/t/tcr39an/tcr39lck3a.dat' u 1:3 with linespoints ls 1 t "TCR",\
'/t/tcr38an/tcr38lck3a.dat' u 1:3 with linespoints ls 2 t "TCR-GOF"
quit
eof
###########################################################
# Use R to get Lck to membrane mean and SD
# Need to accomdate possible 82 columns and fill blanks with NA.
read.table("/t/tcr38an/tcr38lck3b.dat",header=FALSE,col.names=seq(1,39),fill=TRUE)->a38
# Treat data as numeric and strip NA before computation.
nRows=length(a38[,1])
out<-matrix(,nrow=nRows,ncol=3) 
for(i in seq(1,nRows,1)) {
  (i+3)*5-> out[i,1]
  round(mean(as.numeric(a38[i,]),na.rm=TRUE),3)-> out[i,2]
  round((sd(as.numeric(a38[i,]),na.rm=TRUE)/2),3)-> out[i,3]
}
write.table(out,file="tcr38.lck3d.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
# Need to accomdate possible 82 columns and fill blanks with NA.
read.table("/t/tcr39an/tcr39lck3b.dat",header=FALSE,col.names=seq(1,39),fill=TRUE)->a39
# Treat data as numeric and strip NA before computation.
nRows=length(a39[,1])
d<- numeric()
e<- numeric()
out<-matrix(,nrow=nRows,ncol=3) 
for(i in seq(1,nRows,1)) {
  (i+3)*5-> out[i,1]
  round(mean(as.numeric(a39[i,]),na.rm=TRUE),3)-> out[i,2]
  round((sd(as.numeric(a39[i,]),na.rm=TRUE)/2),3)-> out[i,3]
  m<- mean(as.numeric(a38[i,]),na.rm=TRUE ) - mean(as.numeric(a39[i,]),na.rm=TRUE)
  round(m,3)-> d[i]  
  t.test(a38[i,],a39[i,])$p.value-> e[i]
}
write.table(out,file="tcr39.lck3d.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
summary(d)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 -2.133   2.247   5.434   5.244   7.411  13.675
# ^ Median a38-a39 difference = 5.434
# 30/82 = 37% of the frames have a p<0.05. 
###########################################################
# Plot tcr39b and tcr38b Lck to membrane averages 
# Fig 15B
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'tcr.lck4.jpg'
# layout: rows, columns
set multiplot layout 1,2 
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 4
set xtics out nomirror 0,100,430 add ("" 50,"" 150,"" 250,"" 350) offset 0,0.2
set ytics out nomirror 0,5,25 offset 1,0
set xrange [0:430]
set yrange [0:30]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
set ylabel "{/:Bold Lck to membrane (\305)}" offset 2.7,0
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
# tcr39.lck3d
set title "{/:Bold TCR}" offset 0,-1
plot 'tcr39.lck3d.dat' u 1:2:3 with yerrorbars ls 1 t "",\
'' u 1:2 with lines ls 2 t ""
# tcr39.lck3d
set title "{/:Bold TCR-GOF}" offset 0,-1
plot 'tcr38.lck3d.dat' u 1:2:3 with yerrorbars ls 1 t "",\
'' u 1:2 with lines ls 2 t ""
quit
eof
###########################################################
# More extensive Figure 15 C-G plots in tcr38lck7.tcl 

