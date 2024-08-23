# Writes a batch submission shell script to run IQTREE (Nguyen et al. 2014) on a
# linux server running Oracle Grid Engine/Sun Grid Engine.
#
# Arguments: 1) path to .fasta alignment for which you wish to use to run the analysis
#
# Copyright (c) 2024 James Baxter

# import modules
import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(description="a script to write sun grid batch job submission for IQTREE")
    parser.add_argument("file_path")
    args = parser.parse_args()
    return args


def get_filename(input_string):
    match = re.search(r'[^/]+$', input_string)
    if match:
        return match.group(0)
    else:
        return input_string


def main():
    inputs = parse_args()
    relative_path = inputs.file_path
    shell_filename = re.sub('.fasta$', '.sh', relative_path)
    fasta_filename = get_filename(relative_path)



    with open(shell_filename, "w") as file:
        file.write('#!/bin/sh \n')
        file.write('################################################################################\n')
        file.write('################################################################################\n')
        file.write('# This script runs a simple IQTREE run on EDDIE \n')
        file.write('# \n')
        file.write('# It requires the .fasta to be present in the current working directory and \n')
        file.write('# for all trees to require the same command line options.\n')
        file.write('#\n')
        file.write('# Works alongside arrayjobsubmit.sh as part of array job submission, but can\n')
        file.write('# easily be appropriated for submission of a single IQTREE run\n')
        file.write('################################################################################\n')
        file.write('# Grid Engine options\n')
        file.write('#$ -N ' + fasta_filename + '\n')
        file.write('#$ -cwd\n')
        file.write('#$ -pe sharedmem 4\n')
        file.write('#$ -l h_vmem=16G\n')
        file.write('#$ -l h_rt=48:00:00\n')
        file.write('#$ -M james.baxter@ed.ac.uk\n')
        file.write('#$ -P roslin_eeid_aiv\n')
        file.write('#$ -m baes\n')
        file.write('. /etc/profile.d/modules.sh\n')
        file.write('################################################################################\n')
        file.write('module load roslin/iqtree/2.0.5\n')
        file.write('################################################################################\n')
        file.write('# Run the program\n')
        file.write("echo '=============================================='\n")
        file.write("echo '** Hello IQTREE user !**'\n")
        file.write('echo "This job is running on $HOSTNAME"\n')
        file.write("echo 'Start IQTREE with AIV data'\n")
        file.write('iqtree2 -s ' + fasta_filename + ' -bb 1000 -nt 4 \n')
        file.write("echo '** Done **'\n")
        file.write('echo "============================================="\n')
        file.write('################################################################################\n')
        file.write('################################################################################\n')
        file.write('# END #\n')
        file.write('################################################################################\n')
        file.write('################################################################################\n')
        file.close()


if __name__ == "__main__":
    main()
