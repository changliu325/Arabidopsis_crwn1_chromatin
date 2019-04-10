

inputFile=$1
#inputFile: mapped HiC reads to be filtered, it shall has 8 columns:
# chr_a; position_a; direction_a; length_a; chr_b; position_b; direction_b; length_b

cutting_table=$2
#cutting_table: a table recording enzyme cutting sites, it shall has 4 columns:
# chr; start; end; enzyme_id;    the last column is not important

lineNumber=`wc -l $inputFile | awk '{print $1}'`

#if the input file has less than 20 million lines, do following work, otherwise go to the "for" loop:

if [ $lineNumber -lt 20000000 ]
	then

	Rscript filtering_HiC_reads_universal.R $inputFile $cutting_table
else

#in each loop, 20 million mapped read pairs are extracted (~500 Mb) and filtered

	start=1
	for i in $(seq 20000000 20000000 $lineNumber) $lineNumber
		do
			echo Extracting lines ${start} to ${i}
			sed -n ${start},${i}p $inputFile > temp_${inputFile}_line_${start}_to_${i}
			Rscript filtering_HiC_reads_universal.R temp_${inputFile}_line_${start}_to_${i} $cutting_table
			rm temp_${inputFile}_line_${start}_to_${i}
		
			start=$[$i+1]
	done
	
fi

#Lastly, combine all possible results files
cat *.result > HiC_filtered_$inputFile
rm *.result

echo Filtering: ${inputFile} -----Completed!

