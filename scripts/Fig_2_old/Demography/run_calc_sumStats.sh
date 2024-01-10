# Format:
#./calc_sumStats_v2 <alignment_ms_format> <nsamples> <nsites> <begin_sample_pop1> <end_sample_pop1> <begin_sample_pop2> <end_sample_pop2> <tag_pop1> <tag_pop2> <switch_populations> <switch_header>

for file in `ls observed_ms`
do
    sites=`grep 'segsites' observed_ms/$file | sed 's/segsites: //g'`

    ./calc_sumStats_v2 observed_ms/$file 68 $sites 0 16 17 33 AME ASI 0 0 > $file.out
    ./calc_sumStats_v2 observed_ms/$file 68 $sites 0 16 34 50 AME EUR 2 0 >> $file.out
    ./calc_sumStats_v2 observed_ms/$file 68 $sites 0 16 51 67 AME IND 2 0 >> $file.out
    ./calc_sumStats_v2 observed_ms/$file 68 $sites 17 33 34 50 ASI EUR 9 0 >> $file.out
    ./calc_sumStats_v2 observed_ms/$file 68 $sites 17 33 51 67 ASI IND 9 0 >> $file.out
    ./calc_sumStats_v2 observed_ms/$file 68 $sites 34 50 51 67 EUR IND 9 0 >> $file.out
    #echo -e "\n" >> $file.out
done
cat *.out > all_obs_stats.out
