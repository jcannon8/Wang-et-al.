# anEdt1.tcl: Cluster analysis of edt ensemble
# Use vmd -dispdev none
cd /home/jcannon/tcr16/edt
# First goal is to write a trajectory with one protein-only model per frame.
# Make a protein-only psf and pdb for cpptraj to use for clustering.
if {[file exists edt.prot.psf]==0} {
package require psfgen
readpsf /home/jcannon/tcr16/edt/edt.1.psf
coordpdb /home/jcannon/tcr16/edt/edt.1.pdb
mol load psf /home/jcannon/tcr16/edt/edt.1.psf pdb /home/jcannon/tcr16/edt/edt.1.pdb
set del [atomselect top "not protein"]
foreach segid [$del get segid] resid [$del get resid] {
  delatom $segid $resid
}
writepsf edt.prot.psf
writepdb edt.prot.pdb
}
# Because the models are in different molecule files, a protein-only trajectory with
# duplicate frames will be made. Load protein coordinates into frame for each model.
set edtEns [mol new /home/jcannon/tcr16/edt.prot.pdb]
molinfo $edtEns get {numframes frame}
set edtXYZ [atomselect $edtEns protein]
for {set i 1} {$i<=200} {incr i} {
set tcr [mol new /home/jcannon/tcr16/edt/edt.$i.psf]
mol addfile /home/jcannon/tcr16/edt/edt.$i.end.rst type netcdf waitfor all
set mod [atomselect $tcr protein]
# Get model protein coordinates and save into trajectory frame.
set xyz [$mod get {x y z}]
$edtXYZ set {x y z} $xyz
# Add new frame to trajectory (initially duplicate of last).
animate dup $edtEns
# Delete atomselection and mol.
$mod delete 
mol delete $tcr
}
molinfo $edtEns get {numframes frame}
# Save trajectory (frames zero-based)
animate write dcd edtEns.dcd beg 0 end 199 waitfor all $edtEns
quit
##############
# testing
vmd
set tcr [mol new /home/jcannon/tcr16/edt.prot.psf]
mol addfile /home/jcannon/tcr16/edtEns.dcd waitfor all
molinfo $tcr get {numframes frame}
#quit
###########################################################
# Cluster analysis
# Determine Amber residue numbers for CT residues using VMD.
vmd -dispdev none
set edt [mol new /home/jcannon/tcr16/edt.prot.pdb]
[atomselect $edt "chain D and name CA"] get resid
# Last 6JXR chain D is Glu129, CT is Thr130-Lys171
[atomselect $edt "chain D and name CA and resid 130 to 171"] get {resname resid}
[atomselect $edt "chain D and name CA and resid 130 to 171"] get {resname residue}
# ^^ Chain D Thr130-Lys171 is Thr32-Lys73 in Amber residue numbering
# Last 6JXR chain E is Arg155, CT is Lys156-Ile207
[atomselect $edt "chain E and name CA and resid 156 to 207"] get {resname resid}
[atomselect $edt "chain E and name CA and resid 156 to 207"] get {resname residue}
# ^^ Chain E Lys156-Ile207 in Lys105-Ile156 in Amber residue numbering
quit
############################################################
# Analyze the ensemble with 5-25 clusters by two algorithms.
# Use cpptraj for clustering CD3ed CTs
# Chain D CT: Thr130-Lys171 is Amber Thr32-Lys73
# Chain E CT: Lys156-Ile207 is Amber Lys105-Ile156
rm -f means.dat hier.dat
for i in {5..25}; do 
cpptraj<<eof>test.out
parm edt.prot.pdb
trajin edtEns.dcd
# Cluster both chain E and D cytoplasmic CA
cluster out test kmeans clusters $i rms :32-73,105-156@CA
go
eof
echo -n "$i " >> means.dat
grep DBI test.out >> means.dat
done
# 
for i in {5..25}; do 
cpptraj<<eof>test.out
parm edt.prot.pdb
trajin edtEns.dcd
cluster out test hieragglo clusters $i rms :32-73,105-156@CA
go
eof
echo -n "$i " >> hier.dat
grep DBI test.out >> hier.dat
done
rm test test.out
###########################################################
# DBI plot of k-means and hierarchical clustering
gnuplot<<eof
set term jpeg font "arial.ttf,18"
set out 'edtCluster1.jpg' 
set title "CD3ed CT cluster analysis" offset 0,-1 
set xlabel "Cluster number" offset 0,0.5
set ylabel "DBI" offset 1,0
set format y "%.2f"
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
plot 'hier.dat' using 1 : 3 with linespoints ls 1 title "hierarchical", \
'means.dat' using 1 : 3 with linespoints ls 2 title "means"
quit
eof
#
cpptraj<<eof
parm edt.prot.pdb
trajin edtEns.dcd
cluster summary edt14.dat kmeans clusters 16 rms :32-73,105-156@CA
go
eof
# Cluster summary of 14 k-means clusters
cat edt14.dat
#Cluster   Frames     Frac  AvgDist    Stdev Centroid AvgCDist
       0       32    0.160   16.750    2.251       84   15.210
       1       28    0.140   16.961    2.332      171   17.970
       2       24    0.120   17.616    2.270       30   17.824
       3       22    0.110   16.705    2.617      115   22.415
       4       21    0.105   16.569    2.069       35   16.478
       5       18    0.090   17.606    2.649       90   17.635 
       6       15    0.075   16.895    2.371       13   18.697
       7       12    0.060   16.746    2.185       41   17.191
       8       11    0.055   17.108    2.061       32   22.706
       9       10    0.050   16.288    2.092      153   19.708
      10        4    0.020   14.167    1.899      198   20.959
      11        1    0.005    0.000    0.000      187   21.871
      12        1    0.005    0.000    0.000      156   22.675
      13        1    0.005    0.000    0.000      200   22.566
# Pie chart of cluster distribution
gnuplot<<eof
set term jpeg font "arial.ttf,18"
set out 'edt14pie.jpg' 
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
set title "edt CD3ed CT 14 k-means clusters"
plot 'edt14.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
# columns 1-2: x and y coordinates of the center of the disk
# column 3: radius of the disk
# column 4-5: begin and end angles of the region
# column 6: color of the region
eof

