#!/bin/bash
# equilTcr12r.sh: Restrained equilibration of models on Rockfish.
#SBATCH -J tcr33       # job name
#SBATCH -o %j.out      # output and error file name (%j expands to jobID)
#SBATCH -p a100        # Queue name: defq or a100
#SBATCH -N 1
#SBATCH --cpus-per-gpu=1
#SBATCH --gres=gpu:1
#SBATCH -x gpu17          # Exclude the V100 node, gpu01, and gpu17
#SBATCH -A mcb140208_gpu
#SBATCH -t 24:00:00       # run time (hh:mm:ss) - 48 hour (max)

# Load environment variables
module load amber/20
# This job script was used for initial restrained equilibration of many models.
#name=tcr12
#name=CD3edt3
#name=tcr30
name=tcr33
#name=mn34b
#name=tcr35b
#name=tcr35c
#name=tcr33b
#name=tcr33s
#name=tcr35s
#name=mn33c
#name=mn34c
#name=tcr33f
#name=tcr35f
declare -a fc=("250.0" "100.0" "50.0" "50.0" "25.0")
for ((i=0; i<=5;i++)); do
((j=$i+1))
istep="step6.${j}_equilibration"
pstep="$name.eq$i"
if [[ $i == 0 ]] ;then pstep="$name.min";fi
if [[ $i -lt 5 ]] ;then
sed -e "s/FC/${fc[$i]}/g" dihe.restraint > ${istep}.rest
fi
#echo "$i ${fc[$i]} $istep $pstep"
pmemd.cuda -O -i $istep.mdin -p $name.top -c ${pstep}.rst \
-ref $name.crd \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf >& /dev/null

done


