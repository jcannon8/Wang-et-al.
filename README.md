# TCR activation repositions Nck and Lck binding sites Supplementary information 
**Build** Files and scripts to build models  
makeTorsions.cpp: Write tleap impose statements based on HGBLE string.  
makeHGBLE.cpp : Output random HGBLE string based on sequence and frequency table.  
HGBLEfreq2: Secondary structure frequencies for IDR's.  
MakeEpsilon.sh: Create random CD3epsilon cytoplamsic domain models.  
clusterEps.sh: Run clustering on eps ensemble.    
makeCDed4.tcl: Make CD3 epsilon-delta (edt) dimers with random CTs fused to truncated subunits.  
thread2.tcl: Detect threading through Pro, Tyr, His, and Trp sidechains in edt models.  
makeEdn2.tcl: Make CD3 epsilon-delta (edn) dimers with random CTs fused to truncated subunits with Nck bound to CD3 epsilon PRS.  
makeEdn3.tcl: Make CD3 epsilon-delta (edm, edr) dimers with random CTs fused to truncated subunits with Nck bound to CD3 epsilon PRS.  
CD3edt4.psf, CD3edt4.ref2.pdb: Parent of edt dimers used by makeCDed4.tcl.  
cd3edn.psf, cd3edn.ref1.pdb: Parent of edn dimers used by makeEdn2.tcl.   
step5_input.psf, cd3edm.ref1.pdb: Parent of edm, edr dimers used by makeEdn3.tcl. Edm used 2jxb3.c0 Nck-PRS.   
anNck3.sh: Build edr Nck-PRS interaction restraints.

**Run** Job scripts to run MD of models  
tcr33minR.sh: Minimization of tcr33b,c,s and tcr35s on Rockfish.  
edt5min.sh: Minimization of edt.* on Rockfish  
edtNPT4.sh: Heat and two 500 ps NPT of edt models on Expanse  
edtNPT2.sh: continue 5 ns NPT of edt models on Rockfish    
edrNPT4.sh: Heat and two 500 ps NPT of edr models on Expanse   

**Examples** Eight TCR and TCR-GOF structures at 425 ns from figure 14.

**Analyze** Analysis scripts  
anEdt1.tcl: Cluster analysis of edt ensemble
