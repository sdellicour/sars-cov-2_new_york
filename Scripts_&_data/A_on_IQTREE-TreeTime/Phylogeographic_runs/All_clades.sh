#!/bin/bash
#SBATCH --job-name=SARS-CoV-2_NewYork
#SBATCH --time=235:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=10240
#SBATCH --qos=gpu_prio

module load beagle-lib/3.0.2-fosscuda-2018b

cd
cd NCOV
java -jar beast_1104.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite All_clades.xml
