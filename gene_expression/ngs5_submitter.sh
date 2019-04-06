#!/bin/bash
#SBATCH --job-name=job_submitter
#SBATCH --nodes=1
#SBATCH --mem=200MB
#SBATCH --time=1:00:00
#SBATCH --error=job_sub_error.txt
#SBATCH --output=job_sub_stdout.txt

sbatch --array=1 ngs5_script.sh  EV1_R1 EV1_R2
sbatch --array=1 ngs5_script.sh  EV2_R1 EV2_R2
sbatch --array=1 ngs5_script.sh  EV3_R1 EV3_R2
sbatch --array=1 ngs5_script.sh  EV4_R1 EV4_R2
sbatch --array=1 ngs5_script.sh  EV5_R1 EV5_R2
sbatch --array=1 ngs5_script.sh  EV6_R1 EV6_R2
sbatch --array=1 ngs5_script.sh  mVLT1_R1 mVLT1_R2
sbatch --array=1 ngs5_script.sh  mVLT2_R1 mVLT2_R2
sbatch --array=1 ngs5_script.sh  mVLT3_R1 mVLT3_R2
