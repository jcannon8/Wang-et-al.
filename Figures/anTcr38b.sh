# anTcr38b.sh: Analyze KR CT membrane association in trc38a MD (S10 fig).
# vmd -dispdev none on Eire. 
# Remember VMD resid, residue, and index are zero-based. 
# Amber residue and index are one-based.
# tcr38 has chains in ABDEFGMN order
set mod tcr38.0
set tcr [mol new $mod.psf]
mol addfile $mod.eq6.rst type netcdf waitfor all
# Look for Arg and Lys landing on membrane.
# Chain F CT resid 157-207 
set fRK [atomselect $tcr "chain F and resid 157 to 207 and resname ARG LYS and name CA"]
$fRK get {resname residue}; # 13 residues
$fRK get {resname resid}; # 13 residues
# Chain E CT resid 156-207 
set eRK [atomselect $tcr "chain E and resid 156 to 207 and resname ARG LYS and name CA"]
$eRK get {resname residue}; # 13 residues
$eRK get {resname resid}; # 13 residues
# Chain D CT resid 130-171
set dRK [atomselect $tcr "chain D and resid 130 to 171 and resname ARG LYS and name CA"]
$dRK get {resname residue}; # 6 residues
$dRK get {resname resid}; # 6 residues
# Chain G CT resid 139-182
set gRK [atomselect $tcr "chain G and resid 139 to 182 and resname ARG LYS and name CA"]
$gRK get {resname residue}; # 6 residues
$gRK get {resname resid}; # 6 residues
# Get the first and last membrane residues
set mem [atomselect $tcr "chain J"] 
set memI [$mem get residue]
lindex $memI 0; lindex $memI end; # Residues 469-2068
###########################################################
# Analyze Arg, Lys distance to membrane in tcr38 MD
# Run ccptraj to get distance to membrane
mem=":460-2068"
# Important to use "noimage" otherwise it will get the shortest to membrane distance,
# which might be in another periodic box.
ncParm="noimage mindist"
# Modified model list with models with CT that stay within box.
for mod in 0 2 3 5 7 8 9 10 11; do
name=tcr38a.$mod
echo $name
# Build collection of trajectory *.rst inputs
rm -f trajx
echo time>time.dat
# Up to 123 ns (eq28)
# Up to ~200 ns (eq56)
for ((j=3; j<=56;j++)); do
if [ -e /t/tcr38/tcr38.$mod.eq$j.rst ]; then
echo trajin /t/tcr38/tcr38.$mod.eq$j.rst >> trajx
# Get restart file time (ps) from last frame.
grep TIME /t/tcr38/tcr38.$mod.eq$j.out|tail -n1|awk '{printf("%7d\n",$6)}' >>time.dat
fi
done
cpptraj<<eof
# Here the Charmm format topology is being read without problem.
parm /t/tcr38/tcr38.$mod.top
# Load trajectories like "source"
readinput trajx
#
# Minimum distance between Chain F RK and membrane.
# name is 6JXK residue name
nativecontacts name fR179 :282 $mem $ncParm
nativecontacts name fR192 :295 $mem $ncParm
nativecontacts name fK193 :296 $mem $ncParm
nativecontacts name fR196 :299 $mem $ncParm
nativecontacts name fR205 :308 $mem $ncParm
nativecontacts name fR206 :309 $mem $ncParm
# Minimum distance between Chain E RK and membrane.
nativecontacts name eR179 :197 $mem $ncParm
nativecontacts name eR192 :210 $mem $ncParm
nativecontacts name eK193 :211 $mem $ncParm
nativecontacts name eR196 :214 $mem $ncParm
nativecontacts name eR205 :223 $mem $ncParm
nativecontacts name eR206 :224 $mem $ncParm
# Minimum distance between Chain D RK and membrane.
nativecontacts name dR132 :103 $mem $ncParm
nativecontacts name dR144 :115 $mem $ncParm
nativecontacts name dR153 :124 $mem $ncParm
nativecontacts name dK155 :126 $mem $ncParm
nativecontacts name dR169 :140 $mem $ncParm
nativecontacts name dK171 :142 $mem $ncParm
# Minimum distance between Chain G RK and membrane.
nativecontacts name gR143 :346 $mem $ncParm
nativecontacts name gR146 :349 $mem $ncParm
nativecontacts name gK150 :353 $mem $ncParm
nativecontacts name gK164 :367 $mem $ncParm
nativecontacts name gR166 :369 $mem $ncParm
nativecontacts name gR180 :383 $mem $ncParm
nativecontacts name gR181 :384 $mem $ncParm
# This "go" is necessary before the write commands.
go
# Chain F
write memF.dat fR179[mindist] fR192[mindist] fK193[mindist] \
fR196[mindist] fR205[mindist] fR206[mindist]
# Chain E
write memE.dat eR179[mindist] eR192[mindist] eK193[mindist] \
eR196[mindist] eR205[mindist] eR206[mindist]
# Chain D
write memD.dat dR132[mindist] dR144[mindist] dR153[mindist] \
dK155[mindist] dR169[mindist] dK171[mindist]
# Chain G
write memG.dat gR143[mindist] gR146[mindist] gK150[mindist] \
gK164[mindist] gR166[mindist] gR180[mindist] gR181[mindist]
# Just most C-terminal RK
write rk.dat eR206[mindist] fR206[mindist] dK171[mindist] gR180[mindist]
go
eof
####################################
# Use paste to concatenate data to time
# "2" versions used noimage
paste time.dat memF.dat >$name.memF2.dat
paste time.dat memE.dat >$name.memE2.dat
paste time.dat memD.dat >$name.memD2.dat
paste time.dat memG.dat >$name.memG2.dat
paste time.dat rk.dat >$name.mem.dat
rm time.dat memF.dat memE.dat memD.dat memG.dat
rm -f trajx
done
###########################################################
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
# For 123 ns plots, b has imaged distances, c has noimage
#set out 'tcr38.mem123c.jpg'
# For ~200 ns plots
set out 'tcr38.mem200.jpg'
# layout: rows, columns
set multiplot layout 10,4 
set border 3
set rmargin 0
set lmargin 5
set tmargin 0
set bmargin 0
unset xtics
#set ytics out nomirror 0,25,100 offset 1,0 
set ytics out nomirror 0,10,30 offset 1,0 
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
#set xrange [0:140]
set xrange [0:250]
set yrange [0:50]
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "orange" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 6 lc rgb "brown" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 7 lc rgb "cyan" lt 1 lw 3 pt 7 pi -1 ps 1.0
#do for [n = 0:8] {
# Modified model list
do for [n in "0 2 3 5 7 8 9 10 11"] {
#  set title sprintf("tcr19.%s",n) offset 0,-1
#  if (n==147) {set ylabel "{/:Bold Minimum distance to membrane (\305)}" offset 2.5,0}
#  else {set ylabel " "}
  if (n==11) {
    set xlabel "{/:Bold Time (ns)}"
#    set xtics out nomirror 0,40,130
    set xtics out nomirror 0,50,220
  }
  else {unset xlabel}
  # Plot chain E data, label with 6JXR resid
  filename = sprintf("tcr38a.%s.memE2.dat",n)
  plot filename using (\$1/1000) : 3 with lines ls 1 t "R179",\
  '' using (\$1/1000) : 4 with lines ls 2 t "R192",\
  '' using (\$1/1000) : 5 with lines ls 3 t "K193",\
  '' using (\$1/1000) : 6 with lines ls 4 t "R196",\
  '' using (\$1/1000) : 7 with lines ls 5 t "R205",\
  '' using (\$1/1000) : 8 with lines ls 6 t "R206" 
  # Plot chain F data, label with 6JXR resid  
  filename = sprintf("tcr38a.%s.memF2.dat",n)
  plot filename using (\$1/1000) : 3 with lines ls 1 t "R179",\
  '' using (\$1/1000) : 4 with lines ls 2 t "R192",\
  '' using (\$1/1000) : 5 with lines ls 3 t "K193",\
  '' using (\$1/1000) : 6 with lines ls 4 t "R196",\
  '' using (\$1/1000) : 7 with lines ls 5 t "R205",\
  '' using (\$1/1000) : 8 with lines ls 6 t "R206" 
  # Plot chain D data, label with 6JXR resid  
  filename = sprintf("tcr38a.%s.memD2.dat",n)
  plot filename using (\$1/1000) : 3 with lines ls 1 t "R132",\
  '' using (\$1/1000) : 4 with lines ls 2 t "R144",\
  '' using (\$1/1000) : 5 with lines ls 3 t "R153",\
  '' using (\$1/1000) : 6 with lines ls 4 t "R155",\
  '' using (\$1/1000) : 7 with lines ls 5 t "R169",\
  '' using (\$1/1000) : 8 with lines ls 6 t "K171" 
  # Plot chain G data, label with 6JXR resid  
  filename = sprintf("tcr38a.%s.memG2.dat",n)
  plot filename using (\$1/1000) : 3 with lines ls 1 t "R143",\
  '' using (\$1/1000) : 4 with lines ls 2 t "R146",\
  '' using (\$1/1000) : 5 with lines ls 3 t "K150",\
  '' using (\$1/1000) : 6 with lines ls 4 t "K164",\
  '' using (\$1/1000) : 7 with lines ls 5 t "R166",\
  '' using (\$1/1000) : 9 with lines ls 7 t "R181"
}
quit
eof
###########################################################
# Use R to get mean and SD distance for all models at each time
# Files tcr39.*.**.mem.dat have 4 RK to membrane distances in columns 3-7
files <- Sys.glob("tcr38a.*.mem.dat")
nMod=length(files)
# Number of residues
nRes=4
# files have 44-52 frames
nFrame=44
# D array has dimensions model, residue, frame
D<-array(,dim=c(nMod,nRes,nFrame))
for(i in seq(1,nMod,1)) {
  read.table(files[i],skip=1)->a
  # residues are columns starting with column 3.
  for(j in seq(3,(nRes+2),1)) {D[i,(j-2),]<-a[1:nFrame,j]}
}
# Output mean and SD of residue (r) distance for each frame (c) over all models 
M<-matrix(,nrow=nFrame,ncol=nRes)
S<-matrix(,nrow=nFrame,ncol=nRes)
for(c in seq(1,nRes,1)) {
  for(r in seq(1,nFrame,1)) {
    # r=residue index, c=time index
    mean(D[,c,r])->M[r,c]
    sd(D[,c,r])/2->S[r,c]
  }
}
write.table(M,file="tcr38a.memAve.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(S,file="tcr38a.memSD.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)

###########################################################
# combine mean and SD
# Figure S10 fig
cd /t/tcr38an
paste tcr38a.memAve.dat tcr38a.memSD.dat > tcr38a.memAS.dat
# Plot average distances for CD3e 
residues="R206e R206f K171d R181g"
gnuplot<<eof
#set term jpeg font "arial.ttf,18" size 1280,960
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr38a.mem07.jpg'
# layout: rows, columns
#set multiplot layout 2,4 title "{/:Bold TCR^Q, TCR-GOF^Q CD3 CT membrane association}"
set multiplot layout 2,4
set border 3
set rmargin 1
set lmargin 3
set xtics out nomirror 0,50,200
set ytics out nomirror offset 1,0
set grid ytics
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:220]
set yrange [0:30]
set xtics 0,100,200 add ("" 50,"" 150)
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set xlabel "{/:Bold Time (ns)}" offset 0,0.5
set errorbars small
# Time is row number * 5 ns, 
meanCol = 1
SDcol = 5
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-2
  # This puts one ylabel per set of six plots.
#  if (meanCol==1) {set ylabel "{/:Bold Average distance to membrane (\305)}" offset 3,-10}
#  else {set ylabel " "}
  plot 'tcr39a.memAS.dat' u (\$0*5+5):meanCol:SDcol with yerrorbars ls 1 t "",\
  'tcr39a.memAS.dat' u (\$0*5+5):meanCol with lines ls 2 t "" 
  meanCol = meanCol+1
  SDcol = SDcol+1
}
meanCol = 1
SDcol = 5
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-2
  # This puts one ylabel per set of six plots.
#  if (meanCol==1) {set ylabel "{/:Bold Average distance to membrane (\305)}" offset 3,-10}
#  else {set ylabel " "}
  plot 'tcr38a.memAS.dat' u (\$0*5+5):meanCol:SDcol with yerrorbars ls 1 t "",\
  'tcr38a.memAS.dat' u (\$0*5+5):meanCol with lines ls 2 t "" 
  meanCol = meanCol+1
  SDcol = SDcol+1
}
quit
eof
