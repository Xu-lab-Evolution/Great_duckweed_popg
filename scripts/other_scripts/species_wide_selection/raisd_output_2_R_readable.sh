sub="//"
while read -r line;
do
	if [[ "$line" == *"$sub"* ]]
		then
			IFS=' ' read -ra substrings <<< "$line";
			chrom=${substrings[1]};
		else
			echo -e "$chrom\t$line";
	fi
done < $1 #Here goes the RAiSD output Report file
