# Parallel BayesR 

#________________________________________________________
# Parallelise Step 1 
cd /data_dir/ewas/bayesr/gmrm-omi/build_gcc_openmpi
micromamba run -n BayesR mycommand

#Run this first
export OMP_NUM_THREADS=1

#Parallelising function
function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do #counts no background jobs, compares to no. specified when function called
wait -n #waits for one to finish if no. is >= no. specified
done
}


#Run
for i in /home_dir/gmrm/input/*.phen; 
do
A=$(basename -- "$i" .phen)
echo $i
echo $A

./gmrm \
--bin-file /home_dir/gmrm/input/full_scaled_methylation_2602.bin \
--dim-file /home_dir/gmrm/full_dim.dim \
--phen-files $i \
--group-index-file /home_dir/gmrm/full_gri.gri \
--group-mixture-file /home_dir/gmrm/test/test_grm.grm \
--shuffle-markers 1 \ 
--seed 171014 \
--iterations 2000 \
--out-dir /home_dir/gmrm/output/${A}/ &
pwait 3
done
wait
echo "All done" 



#Step 2________________________________________________________________________
cd /data_dir/ewas/bayesr/gmrm-omi/build_gcc_openmpi
micromamba run -n new_BayesR mycommand

#Run this first
export OMP_NUM_THREADS=1

function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do
wait -n
done
}

cd /home_dir/gmrm/input/

for i in *.phen; 
do 
A=$(basename -- "$i" .phen)
echo $i
echo $A

./gmrm \
--bin-file /home_dir/gmrm/input/full_scaled_methylation_2602.bin \
--dim-file /home_dir/gmrm/full_dim.dim \
--phen-files $i \
--iterations 2000 \
--burn-in 750 \
--model linear \
--test \
--in-name-base ../$A \
--out-dir /home_dir/gmrm/output/${A}/new &

pwait 2
done
wait
echo "All done" 