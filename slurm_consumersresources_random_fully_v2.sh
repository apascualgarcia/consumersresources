#!/bin/bash
#SBATCH -J myjob
#SBATCH -o output_file.output
#SBATCH -e error.errors
#SBATCH --mem=5000
#SBATCH --time=10-1:1:1
#SBATCH --cpus-per-task=2

# Example with a job called myjob, with the output redirected to output_file.output
# 50000 Mb to be used, 7 days 1 hour 1 minute and 1 second as maximum execution time.
# python plasticmix_simulation.py > output_mixplastic_second.txt # Execution line, chage it!
#./run_commands commands/study_systems_NR25_NS25_Alberto_2finish.txt 20
#./run_commands commands/study_systems_NR25_NS25_Alberto_optimal_matrix.txt 10
./run_commands commands/study_systems_NR25_NS25_alg-L-BFGS-B_random_fully_v2.txt 2
