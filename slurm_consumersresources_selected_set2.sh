#!/bin/bash
#SBATCH -J crm_sel_set2
#SBATCH -o output_file.output
#SBATCH -e error.errors
#SBATCH --mem=25000
#SBATCH --time=13-1:1:1
#SBATCH --cpus-per-task=30

# Example with a job called myjob, with the output redirected to output_file.output
# 50000 Mb to be used, 7 days 1 hour 1 minute and 1 second as maximum execution time.
# python plasticmix_simulation.py > output_mixplastic_second.txt # Execution line, chage it!
./run_commands commands/study_systems_selected_set2.txt 30
