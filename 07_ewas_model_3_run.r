#Osca run fast-linear, GS full (with protein data). Basic model

#Screen

#--------------------------------------------------------------------------------
#Batch 1
cd /data_dir/ewas/input_files/protein_data/batches/batch1 \

for i in *.phen 

do

  A=$(basename -- "$i" .phen)

echo $i
echo $A


#edit this
time osca_Linux \
--linear \
--befile /data_dir/ewas/input_files/osca_full \
--pheno /data_dir/ewas/input_files/protein_data/batches/batch1/$i \
--qcovar /data_dir/ewas/input_files/protein_data/batches/model/quant_smok.cov \
--fast-linear \
--out /data_dir/ewas/output_files/model_3/fast_linear/${A}_smok \

done
#--------------------------------------------------------------------------------
#Batch2
cd /data_dir/ewas/input_files/protein_data/batches/batch2 \

for i in *.phen 

do

  A=$(basename -- "$i" .phen)

echo $i
echo $A


time osca_Linux \
--linear \
--befile /data_dir/ewas/input_files/osca_full \
--pheno /data_dir/ewas/input_files/protein_data/batches/batch2/$i \
--qcovar /data_dir/ewas/input_files/protein_data/batches/model/quant_smok.cov \
--fast-linear \
--out /data_dir/ewas/output_files/model_3/fast_linear/${A}_smok \

done

#--------------------------------------------------------------------------------
#Batch 3
cd /data_dir/ewas/input_files/protein_data/batches/batch3 \


for i in *.phen 

do

  A=$(basename -- "$i" .phen)

echo $i
echo $A


time osca_Linux \
--linear \
--befile /data_dir/ewas/input_files/osca_full \
--pheno /data_dir/ewas/input_files/protein_data/batches/batch3/$i \
--qcovar /data_dir/ewas/input_files/protein_data/batches/model/quant_smok.cov \
--fast-linear \
--out /data_dir/ewas/output_files/model_3/fast_linear/${A}_smok

done

#--------------------------------------------------------------------------------
#Batch 4

cd /data_dir/ewas/input_files/protein_data/batches/batch4 \


for i in *.phen 

do

  A=$(basename -- "$i" .phen)

echo $i
echo $A


time osca_Linux \
--linear \
--befile /data_dir/ewas/input_files/osca_full \
--pheno /data_dir/ewas/input_files/protein_data/batches/batch4/$i \
--qcovar /data_dir/ewas/input_files/protein_data/batches/model/quant_smok.cov \
--fast-linear \
--out /data_dir/ewas/output_files/model_3/fast_linear/${A}_smok

done