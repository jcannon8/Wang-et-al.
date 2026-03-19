#!/bin/bash                                              
# clusterEps.sh: Run clustering on eps ensemble
cd /home/jcannon/tcr2/cd3f
mv temp.top eps.top
# ff14IDPSFF force field topologies are not understood by cpptraj.
topology=eps.top
# First concatenate last, 5 ns frame into ensemble
# makeEpsilon.sh saved the restart after 5 ns MD.
# Since cpptraj cannot append to a nonexistant file, add the first frame
rm -f epsEns.cdf
cpptraj<<eof
parm $topology
trajin traj/eps.1.rst
# Save one frame in epsEns.cdf for appending below.
trajout epsEns.cdf netcdf
go
eof
# Go through all the trajectory files in cd3d/traj
for f in traj/eps.*.rst; do
# The order is not numerical
echo -n "$f "
cpptraj<<eof
parm $topology
trajin $f $frame $frame
# epsEns.cdf must previously exist, otherwise this fails.
trajout epsEns.cdf netcdf append
go
eof
done
#################################################
# Analyze the ensemble with 5-25 clusters by two algorithms.
# ptraj is deprecated, use cpptraj for clustering
# This compares the whole structure even though PRS is restained.
# The second clustering omits residues PRS residues 30-38 from comparison.
rm -f means.dat hier.dat
for i in {5..25}; do 
cpptraj<<eof>test.out
parm $topology
trajin epsEns.cdf 2 last
cluster out test kmeans clusters $i rms @CA
#cluster out test kmeans clusters $i rms :1-29,39-57@CA
#cluster out test kmeans clusters $i rms :1-29@CA
#cluster out test kmeans clusters $i rms :39-57@CA
go
eof
echo -n "$i " >> means.dat
#grep SSR test.out >> means.dat
#grep pSF test.out >> means.dat
grep DBI test.out >> means.dat
done
 
for i in {5..25}; do 
cpptraj<<eof>test.out
parm $topology
trajin epsEns.cdf 2 last
cluster out test hieragglo clusters $i rms @CA
#cluster out test hieragglo clusters $i rms :1-29,39-57@CA
#cluster out test hieragglo clusters $i rms :1-29@CA
#cluster out test hieragglo clusters $i rms :39-57@CA
go
eof
echo -n "$i " >> hier.dat
#grep SSR test.out >> hier.dat
#grep pSF test.out >> hier.dat
grep DBI test.out >> hier.dat
done
rm test test.out
#################################################
# Figure 11A
cd /home/jcannon/tcr2/cd3f
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 500,375
set out 'epsCluster2.jpg'
set border 3 
set title "{/:Bold CD3ε CT cluster analysis}" offset 0,-1 
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/home/jcannon/tcr2/cd3f/hier.dat' using 1:3 with linespoints ls 1 title "hierarchical", \
'/home/jcannon/tcr2/cd3f/means.dat' using 1:3 with linespoints ls 2  title "means"
quit
eof
###########################################################
topology=eps.top
# Examine kmeans whole clusters.
cpptraj<<eof
parm $topology
trajin epsEns.cdf 2 last
cluster summary whole.dat kmeans clusters 14 rms @CA
go
eof
# Examine kmeans PRSOmit clusters.
cpptraj<<eof
parm $topology
trajin epsEns.cdf 2 last
cluster summary meansOmit.dat kmeans clusters 11 rms :1-29,39-57@CA
go
eof
# Examine kmeans ITAM clusters.
cpptraj<<eof
parm $topology
trajin epsEns.cdf 2 last
cluster summary meansITAM.dat kmeans clusters 13 rms :39-57@CA 
go
eof
# Examine means BRS clusters. 
cpptraj<<eof
parm $topology
trajin epsEns.cdf 2 last
cluster summary meansBRS.dat kmeans clusters 12 rms :1-29@CA 
go
eof
#
cat whole.dat
#Cluster   Frames     Frac  AvgDist    Stdev Centroid AvgCDist
       0       45    0.188   10.401    1.409      102   11.671
       1       38    0.158   11.049    1.549      180   11.180
       2       36    0.150   10.660    1.560       36   14.285
       3       28    0.117   10.718    1.468      120   12.248
       4       25    0.104   10.476    1.410      215   11.593
       5       19    0.079   10.665    1.617       38   11.712
       6       17    0.071   10.821    1.776       22   12.557
       7        7    0.029   10.207    1.135      186   14.593
       8        6    0.025   10.073    1.033       12   14.481
       9        6    0.025   10.069    1.727      223   19.050
      10        6    0.025   11.080    1.218        1   13.531
      11        5    0.021    9.144    1.167      114   13.713
      12        1    0.004    0.000    0.000       46   14.831
      13        1    0.004    0.000    0.000      150   17.835
cat meansOmit.dat
#Cluster   Frames     Frac  AvgDist    Stdev Centroid AvgCDist
       0       49    0.204   10.661    1.545        4   11.036
       1       37    0.154   11.140    1.666      215   12.806
       2       35    0.146   10.869    1.598       36   13.102
       3       32    0.133   10.749    1.486       62   10.178
       4       23    0.096   10.734    1.658      163   11.535
       5       22    0.092   10.443    1.461      120   11.032
       6       15    0.062   10.472    1.582      170   12.935
       7       12    0.050   11.168    1.577       38   11.076
       8       10    0.042   11.043    1.562       82   14.376
       9        3    0.013    9.589    0.756      123   18.621
      10        2    0.008    7.526    0.000       54   12.222
cat meansITAM.dat
#Cluster   Frames     Frac  AvgDist    Stdev Centroid AvgCDist
       0       37    0.154    4.946    0.793       48    5.195
       1       34    0.142    4.976    0.848      239    5.544
       2       26    0.108    4.375    0.782      146    5.300
       3       24    0.100    4.464    0.812      111    6.057
       4       20    0.083    4.946    0.798       51    6.294
       5       20    0.083    4.994    0.758      221    5.650
       6       18    0.075    5.113    0.849      191    5.839
       7       17    0.071    5.017    0.820       77    5.637
       8       14    0.058    5.264    0.753      117    5.436
       9       14    0.058    5.027    0.874       70    6.023
      10        8    0.033    4.580    0.745      184    6.540
      11        7    0.029    4.551    0.655      142    7.996
      12        1    0.004    0.000    0.000      150    7.473
cat meansBRS.dat
#Cluster   Frames     Frac  AvgDist    Stdev Centroid AvgCDist
       0       37    0.154    7.112    1.026       35    7.021
       1       35    0.146    6.288    1.050       26    6.964
       2       34    0.142    6.088    0.903       50    7.155
       3       34    0.142    6.715    1.005      115    7.263
       4       26    0.108    6.458    0.924       96    7.899
       5       25    0.104    6.817    0.985       22    7.086
       6       17    0.071    7.543    1.081      239    7.570
       7       12    0.050    6.243    0.950        1   10.154
       8       10    0.042    6.594    1.034       78    7.721
       9        5    0.021    6.535    0.821      144   10.446
      10        4    0.017    6.809    0.603      176    9.733
      11        1    0.004    0.000    0.000       58    9.823

# AvgDist=Average distance between members within the cluster
# AvgCDist=Average distance of this cluster to every other cluster.
###########################################################
# Make pie charts of summaries above.
# Modified from 
# https://stackoverflow.com/questions/31896718/generation-of-pie-chart-using-gnuplot
gnuplot<<eof
set term jpeg font "arial.ttf,18"
set out 'eps3pie.jpg' 
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0
radius=1
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1            # equal scale length
set xrange [-radius:radius]  # [-1:-1] 
set yrange [-radius:radius]  # [-1:1]
pos = 0            # init angle
color = 0          # init color
#set title "{/:Bold=32 CD3e 14 k-means}"
set title "CD3e CT 14 k-means clusters"
plot '/home/jcannon/tcr2/cd3f/whole.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
# columns 1-2: x and y coordinates of the center of the disk
# column 3: radius of the disk
# column 4-5: begin and end angles of the region
# column 6: color of the region
eof
#
# Use 14 means with whole CT to make references to start 10 ns MD with runEps1.sh
# for 0-11, omitting singletons.
###########################################################
