To generate the input file for 3P-CLR:

1) Run generate_input_file_part1.sh (run it only once). It generates files with frequency counts of each allele (e.g. countsAME.frq.count, etc.)
2) Run get_geneticPosition.R (This has to be run for each chromosome separately. It will generate an output genetic_positions_chrN.txt where N is the chromosome in question).
3) Run generate_input_file_part2.R (Run this also per chromosome. It generates the final input file for 3P-CLR).
4) Run generate_input_drift.R (Run it also per chromosome. It uses the input from the previous step. It generates an input file for the script CalcDriftsF3.R).
5) Run CalcDriftsF3.R (This script is part of 3P-CLR. Run it per chromosome. It generates the drift parameters to be passed to 3P-CLR. A readme is available for this script).

To run 3P-CLR:
6) In the cluster, run launch_chrXX.sh (where XX is the chromosome to be analyzed). You can run all chromosomes in parallel with scripts like this one). An example input file (generated with steps 1 through 5) is given.

To generate/analyze the output files:
7) After all output files are generated in the cluster, run analyze_output.R to produce the final output files.
8) Run also the bed2genes.sh script to obtain gene IDs from the selected regions produced by 3P-CLR. 
