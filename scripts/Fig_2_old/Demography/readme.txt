Demography analyses were done in several steps.
First, intergenic loci were selected from all chromosomes and their sequences converted to ms-format with fasta2ms.c
Second, summary statistics from each intergenic locus were calculated with calc_sumStats_v2.c and mean and variances calculated with get_mean_var_obs.R
The resulting vector constitutes the "observed" summary statistics.
Then, models were simulated with sim_ms_model_X.sh and summary statistics calculated with calc_sumStats_v2.c and summarized with get_mean_var.R
Subsequently, ABC model choice was performed with the script modelchoice.R
Finally, ABC parameter estimation on the most probable model was performed with the script runABC.R
