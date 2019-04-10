inputFile=$1
# $1: the raw fastq file, either read1 or read2.

bowtieDir=/home-link/bcglc01/bowtie2-2.2.4

#setting of i is flexible, it determines how many bp shall be trimmed from the 3'
#NOTE: some lines only have 4 columns, in such cases, a "*" is added as the 5th column of the temporal file
#scan the file, and determine the mapping (only unique mapping result is considered)
#the final out file has 6 columns. the first column tells you if the read has been uniquely mapped (labeled with "1") or not ("-1")


for i in {75..0..5}
	do
		
		${bowtieDir}/bowtie2 -3 $i --no-hd --no-sq -p 24 --reorder -D 20 -R 5 -N 0 -L 20 -i S,1,0.50 -x /home-link/bcglc01/mapping/index_genome_for_bowtie2/hybrid_ColxLyrata/Ath_index $inputFile -S ${inputFile}_mapping_raw_with_${i}bp_trimmed.sam 
		cut -f 2-4,10,12-13 ${inputFile}_mapping_raw_with_${i}bp_trimmed.sam | awk 'NF==6{print $0;next}{print $0"\t""*""\t""*"}' | awk '($1==16){print $1"\t"$2"\t"length($4)+$3-1"\t"$4"\t"$5"\t"$6;next}{print $0}' | awk '($1==0||$1==16)&&(substr($6,1,2)!="XS"||(substr($5,6)-substr($6,6))>10){print "1\t"$2"\t"$3"\t"$1;next}{print "-1\t"$2"\t"$3"\t"$1}' >${inputFile}_mapping_result_with_${i}bp_trimmed.txt
		rm ${inputFile}_mapping_raw_with_${i}bp_trimmed.sam 
		
	done

# in the mapping result file, each row has four columns: mark; chromosome; position; direction

#combine all mapping results, and look for the best mapping output (the one has the longest enquiry wins)
#if this read is not mapped, results from the first 4 columns are printed instead.

fileList=""
for i in {75..0..5}
	do
	  fileList="${fileList} ${inputFile}_mapping_result_with_${i}bp_trimmed.txt"
	done
	

paste $fileList | awk  'function FindBestMapping() {if ($1==1&&$5==-1){print $1"\t"$2"\t"$3"\t"$4"\t"26;next} else if ($5==1&&$9==-1){print $5"\t"$6"\t"$7"\t"$8"\t"31;next} else if ($9==1&&$13==-1){print $9"\t"$10"\t"$11"\t"$12"\t"36;next} else if ($13==1&&$17==-1){print $13"\t"$14"\t"$15"\t"$16"\t"41;next} else if ($17==1&&$21==-1){print $17"\t"$18"\t"$19"\t"$20"\t"46;next} else if ($21==1&&$25==-1){print $21"\t"$22"\t"$23"\t"$24"\t"51;next} else if ($25==1&&$29==-1){print $25"\t"$26"\t"$27"\t"$28"\t"56;next} else if ($29==1&&$33==-1){print $29"\t"$30"\t"$31"\t"$32"\t"61;next} else if ($33==1&&$37==-1){print $33"\t"$34"\t"$35"\t"$36"\t"66;next} else if ($37==1&&$41==-1){print $37"\t"$38"\t"$39"\t"$40"\t"71;next} else if ($41==1&&$45==-1){print $41"\t"$42"\t"$43"\t"$44"\t"76;next} else if ($45==1&&$49==-1){print $45"\t"$46"\t"$47"\t"$48"\t"81;next} else if ($49==1&&$53==-1){print $49"\t"$50"\t"$51"\t"$52"\t"86;next} else if ($53==1&&$57==-1){print $53"\t"$54"\t"$55"\t"$56"\t"91;next} else if ($57==1&&$61==-1){print $57"\t"$58"\t"$59"\t"$60"\t"96;next} else if ($61==1){print $61"\t"$62"\t"$63"\t"$64"\t"101;next} else print $1"\t"$2"\t"$3"\t"$4"\t"26}{FindBestMapping()}' >  ${inputFile}_mapped

#this file has to be merged with its partner (i.e., the other read)
#the fifth column means sequence length with which mapping result was obtained.
#lastly, remove all temporal mapping files.
rm -f $fileList
