# anTcr39itam1.sh: Compare ITAM Tyr to membrane distances for tcr39*.** and tcr38*.** (Fig 13E,F).
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# https://cran.r-project.org/doc/manuals/r-release/R-intro.html
#
# Files have distances in rows, ITAM Tyr in columns.
cd /home/wangly/tyr_analysis
# tcr38b data files from anTcr38Y2.sh:
# ITAMTYRtoMEM_TCR38eq45.dat
# ITAMTYRtoMEM_TCR38eq85.dat
#
# tcr39b data files from allY_selection2.sh:
# ITAMTYRtoMEM_NoImage45.dat
# ITAMTYRtoMEM_NoImage85.dat
###########################################################
# Use R to get comparative statistics 
# 205 ns comparison using eq45 frame
read.table("/home/wangly/tyr_analysis/NoImage_Test/ITAMTYRtoMEM_NoImage45.dat", header=TRUE)->t39
read.table("/home/wangly/tyr_analysis/TCR38_TYRAnlyz/tcr38_eq45/ITAMTYRtoMEM_TCR38eq45.dat",header=TRUE)->t38
#
p <- numeric()
m38 <- numeric()
m39 <- numeric()
for(i in seq(1,20,1)) {
  # two-sided t-test
  t.test(t39[,i],t38[,i])->q
  p[i] <- format(q$p.value, digits =3, scientific = TRUE)
  m39[i] <- round(mean(t39[, i]),digits=3) 
  m38[i] <- round(mean(t38[, i]),digits=3)
}
n<-c(1:20)
# Output: <ITAM #><Tcr38b mean><Tcr39b mean><p-value>
d<-data.frame(n,m38, m39, p)
write.table(d,"tcrTyr45.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
# 405 ns comparison using eq85 frame
read.table("/home/wangly/tyr_analysis/NoImage_Test/ITAMTYRtoMEM_NoImage85.dat", header=TRUE)->t39
read.table("/home/wangly/tyr_analysis/TCR38_TYRAnlyz/tcr38_eq45/ITAMTYRtoMEM_TCR38eq85.dat",header=TRUE)->t38
#
p <- numeric()
m38 <- numeric()
m39 <- numeric()
for(i in seq(1,20,1)) {
  # two-sided t-test
  t.test(t39[,i],t38[,i])->q
  p[i] <- format(q$p.value, digits =3, scientific = TRUE)
  m39[i] <- round(mean(t39[, i]),digits=3) 
  m38[i] <- round(mean(t38[, i]),digits=3)
}
# Output: <ITAM #><Tcr38b mean><Tcr39b mean><p-value>
d<-data.frame(n,m38, m39, p)
write.table(d,"tcrTyr85.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
q()
###########################################################
# Plot with linear Y-axis with log(p-value)
# Fig 13E, F
cd /home/wangly/tyr_analysis
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 2000,750
set output 'tcr39.tyrSig.jpg'
# layout: rows, columns
set multiplot layout 1,2
set border 3
set rmargin 1
set lmargin 3
#set ylabel 'log(p-value)'
set xtics rotate by -45 nomirror
set ytics -6,1,0 nomirror
set grid ytics
set xrange [0.5:21]
set yrange [-6:0]
set key bottom left center textcolor variable samplen -1 font "arial.ttf,18"
# These have a variety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 2.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 2.0
# Show significance threshold
set arrow from first 0.5,log10(0.05) to first 21,log10(0.05) nohead lt 2
set style line 5 lc rgb "black" lt 1 lw 2 pt 7 pi -1 ps 1.0
set xtics ('AY72' 1, 'AY83' 2, 'AY111' 3, 'AY123' 4, 'AY142' 5, 'AY153' 6, \
'BY72' 7, 'BY83' 8, 'BY111' 9, 'BY123' 10, 'BY142' 11, 'BY153' 12, 'DY149' 13,\
'DY160' 14, 'EY188' 15, 'EY199' 16, 'FY188' 17, 'FY199' 18, 'GY160' 19, 'GY171' 20)
# Make these xtics labels a little smaller so they do not fuse.
set xtics font "arial.ttf,18"
# Use blue and red styles to distinguish relative means.
set title "{/:Bold Comparison at 225 ns}" offset 0,-1.5
file="/home/wangly/tyr_analysis/tcrTyr45.dat"
plot file using 1:($2>$3?log10($4):1/0) with points ls 2 title "TCR-GOF > TCR",\
'' using 1:($2<$3?log10($4):1/0) with points ls 1 title "TCR-GOF < TCR"
#
set title "{/:Bold Comparison at 425 ns}" offset 0,-1.5
file="/home/wangly/tyr_analysis/tcrTyr85.dat"
replot
eof
