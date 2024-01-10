counter=1
calc_sumStats="/scratch/tmp/pduchenb/Duckweed/Simulations/calc_sumStats_v2"
msms_exec="/home/p/pduchenb/Programs/msms/bin/msms"
N_msms=10000
mut=0.0000000002356

while read -r line
do
    log_Ne="$(cut -f 1 -d "," <<< $line)"
    N_ASI_prop="$(cut -f 2 -d "," <<< $line)"
    N_EUR_prop="$(cut -f 3 -d "," <<< $line)"
    N_IND_prop="$(cut -f 4 -d "," <<< $line)"
    log_migE="$(cut -f 5 -d "," <<< $line)"
    log_migI="$(cut -f 6 -d "," <<< $line)"
    log_migA="$(cut -f 7 -d "," <<< $line)"
    timebn="$(cut -f 8 -d "," <<< $line)"
    timecoal_EUR="$(cut -f 9 -d "," <<< $line)"
    timecoal_IND="$(cut -f 10 -d "," <<< $line)"
    timecoal_ASIA="$(cut -f 11 -d "," <<< $line)"
    SAA="$(cut -f 12 -d "," <<< $line)"
    Shet="$(cut -f 13 -d "," <<< $line)"

    Ne="$(cut -f 14 -d "," <<< $line)"
    migE="$(cut -f 15 -d "," <<< $line)"
    migI="$(cut -f 16 -d "," <<< $line)"
    migA="$(cut -f 17 -d "," <<< $line)"
    
    
    N_migE=$(echo "scale = 6; 4 * $Ne * $migE" | bc)
    N_migI=$(echo "scale = 6; 4 * $Ne * $migI" | bc)
    N_migA=$(echo "scale = 6; 4 * $Ne * $migA" | bc)

    echo "Sim $counter --- $log_Ne, $migE, $migI"

    i=1
    while read -r lineRecomb
    do
	rho="$(cut -f 2 -d "," <<< $lineRecomb)" #This is already 4Nr estimated by ldhat, taking into account the locus length.
	length="$(cut -f 3 -d "," <<< $lineRecomb)"
        theta=$(echo "scale = 6; 4 * $Ne * $mut * $length" | bc)
        $msms_exec -N $N_msms -ms 68 1 -t $theta -r $rho $length -I 4 17 17 17 17 -n 2 $N_ASI_prop -n 3 $N_EUR_prop -n 4 $N_IND_prop -ma x $N_migA $N_migE $N_migI $N_migA x $N_migE $N_migI $N_migE $N_migE x 0 $N_migI $N_migI 0 x -SI $timecoal_EUR 4 0.1 0.1 0.1 0.1 -SAA $SAA -SaA $Shet -en $timebn 4 1 -en $timebn 3 1 -en $timebn 2 1 -en $timebn 1 1 -ej $timecoal_EUR 3 2 -ej $timecoal_IND 4 2 -ej $timecoal_ASIA 1 2 > out.ms.tmp.$1
    
	segsites=`grep 'segsites' out.ms.tmp.$1 | sed 's/segsites: //g'`

        #head -n -1 out.ms.tmp.$1 | tail -n +4 > out.ms.$1 #This line (or the next one) make sure the output of msms is understood by calc_sumStats_v2
	tail -n +4 out.ms.tmp.$1 | sed '$d' > out.ms.$1 #Same as last line but for mac. Make sure there's no empty last line in out.ms.$1
	rm out.ms.tmp.$1
	
	$calc_sumStats out.ms.$1 68 $segsites 0 16 17 33 AME ASI 0 0 > sim.out.$1.$i
	$calc_sumStats out.ms.$1 68 $segsites 0 16 34 50 AME EUR 2 0 >> sim.out.$1.$i
	$calc_sumStats out.ms.$1 68 $segsites 0 16 51 67 AME IND 2 0 >> sim.out.$1.$i
	$calc_sumStats out.ms.$1 68 $segsites 17 33 34 50 ASI EUR 9 0 >> sim.out.$1.$i
	$calc_sumStats out.ms.$1 68 $segsites 17 33 51 67 ASI IND 9 0 >> sim.out.$1.$i
	$calc_sumStats out.ms.$1 68 $segsites 34 50 51 67 EUR IND 9 0 >> sim.out.$1.$i
	echo "" >> sim.out.$1.$i
	i=$((i+1))
    done < recomb_table_intergenic_noHeader.csv
    cat sim.out.$1.* > all_stats.out.$1
    Rscript get_mean_var.R all_stats.out.$1 >> sim.out.meanvar.$1
    echo "" >> sim.out.meanvar.$1
    rm sim.out.$1.* all_stats.out.$1
    counter=$((counter+1))
done < $1

