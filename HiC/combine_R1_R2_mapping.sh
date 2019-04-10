# This script combines mapping read1/read2 results from iterative mapping, and retains read pairs with both reads uniquely mapped:
# This script also sorts mapping results and remove PCR duplicates
# NOTE columns of input file shall in such order:
# 	1) mark, 1 or -1
# 	2) chromosome, 
# 	3) position, 
# 	4) direction,
#   5) this column tells you the length of reads mapped. Either present or absent on both files



#
fileA=$1
fileB=$2
touch log_file_mapping_${fileA}.txt
paste $fileA $fileB | awk '$1==1&&$6==1{print $2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10}' > total_mapped_${fileA}

echo 'total mapped reads:' >> log_file_mapping_${fileA}.txt
wc -l total_mapped_${fileA} >> log_file_mapping_${fileA}.txt

awk '$1>$5{print $5"\t"$6"\t"$7"\t"$8"\t"$1"\t"$2"\t"$3"\t"$4;next}{print $0}' total_mapped_${fileA} | awk '$2>$6{print $5"\t"$6"\t"$7"\t"$8"\t"$1"\t"$2"\t"$3"\t"$4;next}{print $0}' | awk '$3>$7{print $5"\t"$6"\t"$7"\t"$8"\t"$1"\t"$2"\t"$3"\t"$4;next}{print $0}' | sort -k1,1 -k2,2 -k3,3n | uniq > final_mapped_unique_${fileA}.txt

echo 'total mapped reads after removing PCR duplicates:' >> log_file_mapping_${fileA}.txt
wc -l final_mapped_unique_${fileA}.txt >> log_file_mapping_${fileA}.txt

#remove the temporal file:
rm total_mapped_${fileA}


#final file has 8 columns as following:
#
#  chr_a; position_a; direction_a; read_length_a; chr_b; position_b; direction_b; read_length_b
#
#NOTE: for direction, "0" means mapped to + strand; "16" means mapped to - strand
