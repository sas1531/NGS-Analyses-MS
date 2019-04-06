#!/bin/bash
#SBATCH --job-name=job_submitter
#SBATCH --nodes=1
#SBATCH --mem=200MB
#SBATCH --time=1:00:00
#SBATCH --error=job_sub_error.txt
#SBATCH --output=job_sub_stdout.txt

sbatch --array=1 ngs6_script.sh  SRR7992458  iCTRL1
sbatch --array=1 ngs6_script.sh  SRR7992461  iCTRL2
sbatch --array=1 ngs6_script.sh  SRR7992460  iCTRL3
sbatch --array=1 ngs6_script.sh  SRR7992450  iDS1
sbatch --array=1 ngs6_script.sh  SRR7992457  iDS2
sbatch --array=1 ngs6_script.sh  SRR7992456  iDS3
sbatch --array=1 ngs6_script.sh  SRR7992455  mCTRL1
sbatch --array=1 ngs6_script.sh  SRR7992454  mCTRL2
sbatch --array=1 ngs6_script.sh  SRR7992459  mCTRL3
sbatch --array=1 ngs6_script.sh  SRR7992453  mDS1
sbatch --array=1 ngs6_script.sh  SRR7992452  mDS2
sbatch --array=1 ngs6_script.sh  SRR7992451  mDS3
