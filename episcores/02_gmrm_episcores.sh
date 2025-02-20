#episcore run

#Step 1

# Parallelise Step 1
cd /data_dir/ewas/new_bayesr/gmrm-omi/build_gcc_openmpi
micromamba run -n new_BayesR mycommand

#Run this first
export OMP_NUM_THREADS=1


#Parallelising function
function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do #counts no background jobs, compares to no. specified when function called
wait -n #waits for one to finish if no. is >= no. specified
done
}


for i in /home_dir/gmrm_scores/input/*.phen; 
do
A=$(basename -- "$i" .phen)
echo $i
echo $A

./gmrm \
--bin-file /home_dir/gmrm_scores/input/scores_methylation.bin \
--dim-file /home_dir/gmrm_scores/input/scores_dim.dim \
--phen-files $i \
--group-index-file /home_dir/gmrm_scores/input/full_gri.gri \
--group-mixture-file /home_dir/gmrm_scores/input/grm.grm \
--shuffle-markers 1 \ 
--seed 171014 \
--iterations 2000 \
--out-dir /home_dir/gmrm_scores/output/${A}/ &
pwait 3
done
wait
echo "All done" 



#Step 2
cd /data_dir/ewas/new_bayesr/gmrm-omi/build_gcc_openmpi
micromamba run -n new_BayesR mycommand

#Run this first
export OMP_NUM_THREADS=1


#Parallelising function
function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do #counts no background jobs, compares to no. specified when function called
wait -n #waits for one to finish if no. is >= no. specified
done
}


for i in /home_dir/gmrm_scores/input/*.phen; 
do 
A=$(basename -- "$i" .phen)
echo $i
echo $A

./gmrm \
--bin-file /home_dir/gmrm_scores/input/scores_methylation.bin \
--dim-file /home_dir/gmrm_scores/input/scores_dim.dim \
--phen-files $i \
--iterations 2000 \
--burn-in 750 \
--model linear \
--test \
--in-name-base $A \
--out-dir /home_dir/gmrm_scores/output/${A} &

pwait 2
done
wait
echo "All done" 