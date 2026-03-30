# anTcr39m.sh: Analyze KR membrane association CT residues in tcr39b (S11 fig).
# derived from anTcr38m.sh
cd /t/tcr39an
# Use VMD to get 1-based Amber residue numbers for CD3 CT KR residues. 
vmd -dispdev none
set tcr [mol new /t/tcr39b/tcr39.0.00.psf]
# Start with C-terminal RK for chains EFDG
# Chain E CT resid 156-207, chain D resid 130-171
set eRK [atomselect $tcr "chain E and resid 156 to 207 and resname ARG LYS and name CA"]
$eRK get {resname resid residue}
# Chain E R206 is 0-based residue R439, 1-based residue R440.
set fRK [atomselect $tcr "chain F and resid 156 to 207 and resname ARG LYS and name CA"]
$fRK get {resname resid residue}
# Chain F R206 is 0-based residue R524, 1-based residue R525.
set dRK [atomselect $tcr "chain D and resid 130 to 171 and resname ARG LYS and name CA"]
$dRK get {resname resid residue}
# Chain D K171 is 0-based residue K357, 1-based residue K358.
set gRK [atomselect $tcr "chain G and resid 139 to 182 and resname ARG LYS and name CA"]
$gRK get {resname resid residue}
# Chain G R181 is 0-based residue R599, 1-based residue R600.
# Next get chain AB KR
set aRK [atomselect $tcr "chain A and resid 92 104 136 164 and resname ARG LYS and name CA"]
$aRK get {resname resid residue}
# Chain A R92, K104, K136, R164 are 0-based R70, K82, K114, R142; 
# 1-based R71, K83, K115, R143
set bRK [atomselect $tcr "chain B and resid 92 104 136 164 and resname ARG LYS and name CA"]
$bRK get {resname resid residue}
# Chain B R92, K104, K136, R164 are 0-based R211, K223, K255, R283; 
# 1-based R213, K225, K260, R288
#
# Get the first and last membrane residues
set mem [atomselect $tcr "chain J"] 
set memI [$mem get residue]
lindex $memI 0; lindex $memI end; # Residues 685-2284 in Amber 1-based numbers
# These membrane residue numbers are the same for all tcr39b models.
quit
residues="R206e R206f K171d R181g R92a K104a K136a R164a R92b K104b K136b R164b"
###########################################################
# Use cpptraj to get mimimum distance to membrane every 5 ns.
mem=":685-2284"
# Important to use "noimage" otherwise it will get the shortest to 
# membrane distance, which might be in another periodic box.
ncParm="noimage mindist"
#for i in 0 2 3 5 7 8 9 10 11; do
for ((i=0;i<=11;i++)); do
for j in 00 01 02 10 11 12 20 21 22; do
name=tcr39.$i.$j
p=/t/tcr39b/$name.top
if [ -e $p ]; then
# Valid tcr39b model
out=$name.mem.dat
# Build collection of trajectory *.rst inputs
rm -f trajx
echo time>time.dat
# Up to 206 ns (eq45)
# Up to 405 ns (eq85)
for ((k=4; k<=85;k++)); do
# Some trajectories have 81 and some have 82 files.
if [ -e /t/tcr39b/$name.eq$k.rst ]; then
echo trajin /t/tcr39b/$name.eq$k.rst >> trajx
# Get restart file time (ps) from last frame.
grep TIME /t/tcr39b/$name.eq$k.out|tail -n1|awk '{printf("%7d\n",$6)}' >>time.dat
fi
done
#
cpptraj<<eof
parm $p
# Load trajectories like "source"
readinput trajx
# 12 RK-membrane distances
nativecontacts name R206e :440 $mem $ncParm
nativecontacts name R206f :525 $mem $ncParm
nativecontacts name K171d :358 $mem $ncParm
nativecontacts name R181g :600 $mem $ncParm
nativecontacts name R92a :71 $mem $ncParm
nativecontacts name K104a :83 $mem $ncParm
nativecontacts name K136a :115 $mem $ncParm
nativecontacts name R164a :143 $mem $ncParm
nativecontacts name R92b :211 $mem $ncParm
nativecontacts name K104b :223 $mem $ncParm
nativecontacts name K136b :255 $mem $ncParm
nativecontacts name R164b :283 $mem $ncParm
# This "go" is necessary before the write command.
go
write rk.dat R206e[mindist] R206f[mindist] K171d[mindist] R181g[mindist] \
R92a[mindist] K104a[mindist] K136a[mindist] R164a[mindist] \
R92b[mindist] K104b[mindist] K136b[mindist] R164b[mindist]
go
eof
# Paste together time and RK distances.
paste time.dat rk.dat >$out
# Done with this model
fi
# Next model
done
done
rm -f trajx time.dat rk.dat
###########################################################
# Use R to get mean and SD distance for all models at each time
# Files tcr39.*.**.mem.dat have 12 RK to membrane distances in columns 3,4,5,..
files <- Sys.glob("tcr39.*.**.mem.dat")
nMod=length(files)
read.table(files[1],skip=1)->a
nRes=12
#nFrame=length(a[,1])
# Some trajectories have 81 and some have 82 frames.
nFrame=81
# D array has dimensions model, residue, frame
D<-array(,dim=c(nMod,nRes,nFrame))
for(i in seq(1,nMod,1)) {
  read.table(files[i],skip=1)->a
  # residues are columns starting with column 3.
  for(j in seq(3,14,1)) {D[i,(j-2),]<-a[1:nFrame,j]}
}
# Output distributions of last frame
write.table(D[,,nFrame],file="tcr39.mem81.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
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
write.table(M,file="tcr39.memAve.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(S,file="tcr39.memSD.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)

###########################################################
# combine mean and SD
paste tcr39.memAve.dat tcr39.memSD.dat > tcr39.memAS.dat
# Plot average RK distances to membrane, S11 fig A 
residues="R206e R206f K171d R181g R92a K104a K136a R164a R92b K104b K136b R164b"
gnuplot<<eof
#set term jpeg font "arial.ttf,18" size 1280,960
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr39.mem07.jpg'
# layout: rows, columns
set multiplot layout 3,4 title "{/:Bold TCR RK membrane association}"
set border 3
set rmargin 1
set lmargin 3
set xtics out nomirror 0,100,400
set ytics out nomirror offset 1,0
set grid ytics
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:420]
set yrange [0:30]
set xtics 0,200,400 add("" 100,"" 300)
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set xlabel "{/:Bold Time (ns)}" offset 0,0.5
set errorbars small
# Time is row number * 5 ns, 
meanCol = 1
SDcol =13
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-2
  # This puts one ylabel per set of six plots.
#  if (meanCol==1) {set ylabel "{/:Bold Distance to membrane (\305)}" offset 3,-10}
#  else {set ylabel " "}
  plot 'tcr39.memAS.dat' u (\$0*5+5):meanCol:SDcol with yerrorbars ls 1 t "",\
  'tcr39.memAS.dat' u (\$0*5+5):meanCol with lines ls 2 t "" 
  meanCol = meanCol+1
  SDcol = SDcol+1
}
quit
eof
###########################################################
# boxplot graphs RK membrane distance distributions for 425 ns frame
# S11 fig C  
cp /t/tcr38an/tcr38.mem81.dat .
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750 
set out 'tcr39.mem08.jpg'
set title "{/:Bold RK Membrane distances }"
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 8 
set xtic nomirror rotate by -45 offset -1 font "arial.ttf,18"
set ytic out nomirror
set yrange [0:50]
set grid ytics
#set ylabel "{/:Bold Distance (Å)}" offset 1,0
set style boxplot nooutliers pointtype 7
# Not showing outliers in tcr39.mem08.jpg.
set style data boxplot
set boxwidth 0.5
set key left horizontal textcolor variable samplen -1 font "arial.ttf,18" 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set xtics ("R206e" 1.5, "R206f" 4.5, "K171d" 7.5, "R181g" 10.5, "R92a" 13.5, "K104a" 16.5,\
"K136a" 19.5, "R154a" 22.5, "R92b" 25.5, "K104b" 28.5, "K136b" 31.5, "R164b" 34.5)
# 
plot 'tcr39.mem81.dat' u (1):1 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (2):1 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (4):2 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (5):2 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (7):3 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (8):3 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (10):4 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (11):4 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (13):5 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (14):5 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (16):6 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (17):6 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (19):7 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (20):7 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (22):8 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (23):8 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (25):9 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (26):9 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (28):10 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (29):10 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (31):11 with boxplot ls 1 t "",\
'tcr38.mem81.dat' u (32):11 with boxplot ls 2 t "",\
'tcr39.mem81.dat' u (34):12 with boxplot ls 1 t "TCR",\
'tcr38.mem81.dat' u (35):12 with boxplot ls 2 t "TCR-GOF"
# Outliers shown in tcr39.mem09.jpg
set out 'tcr39.mem09.jpg'
set style boxplot outliers pointtype 7
set yrange [0:100]
replot
eof
###########################################################
# Use R to get t-test for 425 ns RKmem distances
read.table("tcr38.mem81.dat")->t38
read.table("tcr39.mem81.dat")->t39
t.test(t39[,9],t38[,9])
p <- numeric()
m38 <- numeric()
m39 <- numeric()
dis <- numeric()
for(i in seq(1,12,1)) {
  # two-sided t-test
  t.test(t39[,i],t38[,i])->q  
  p[i] <- format(q$p.value, digits =3)
  m39[i] <- round(mean(t39[, i]),digits=3) 
  m38[i] <- round(mean(t38[, i]),digits=3)
  dis[i] <- round((mean(t38[, i]) - mean(t39[, i])),digits=3) 
}
n<-c("R206e","R206f","K171d","R181g","R92a","K104a","K136a", "R154a", "R92b", "K104b", "K136b","R164b")
d<-data.frame(n,m39,m38,p,dis)
write.table(d,"tcrRKp.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
cat tcrRKp.dat
R206e 7.564 9.745 0.263 2.181
R206f 9.591 6.133 0.0125 -3.458
K171d 12.311 10.998 0.483 -1.313
R181g 6.648 8.681 0.151 2.033
R92a 7.012 13.314 0.000964 6.301
K104a 9.802 14.389 0.0219 4.587
K136a 10.863 14.067 0.197 3.204
R154a 13.822 15.756 0.416 1.934
R92b 7.419 14.589 4.08e-05 7.17
K104b 8.747 13.157 0.018 4.409
K136b 10.05 14.768 0.0247 4.719
R164b 16.05 14.443 0.366 -1.607
