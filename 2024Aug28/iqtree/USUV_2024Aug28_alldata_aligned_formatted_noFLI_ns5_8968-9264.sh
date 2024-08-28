#!/bin/sh 
################################################################################
################################################################################
# This script runs a simple IQTREE run on EDDIE 
# 
# It requires the .fasta to be present in the current working directory and 
# for all trees to require the same command line options.
#
# Works alongside arrayjobsubmit.sh as part of array job submission, but can
# easily be appropriated for submission of a single IQTREE run
################################################################################
# Grid Engine options
#$ -N USUV_2024Aug28_alldata_aligned_formatted_noFLI_ns5_8968-9264.fasta
#$ -cwd
#$ -pe sharedmem 4
#$ -l h_vmem=16G
#$ -l h_rt=48:00:00
#$ -l rl9=FALSE
#$ -M james.baxter@ed.ac.uk
#$ -P roslin_eeid_aiv
#$ -m baes
. /etc/profile.d/modules.sh
################################################################################
module load roslin/iqtree/2.0.5
################################################################################
# Run the program
echo '=============================================='
echo '** Hello IQTREE user !**'
echo "This job is running on $HOSTNAME"
echo 'Start IQTREE with AIV data'
iqtree2 -s USUV_2024Aug28_alldata_aligned_formatted_noFLI_ns5_8968-9264.fasta -bb 1000 -nt 4 -czb
echo '** Done **'
echo "============================================="
################################################################################
################################################################################
# END #
################################################################################
################################################################################
