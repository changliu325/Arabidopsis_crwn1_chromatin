file_R1=$1
#extract the file name (in case directory is included in input)
inputFileName=${file_R1##*/}
#input files, the raw fastq files of read1 

module load devel/python/2.7.12

#mapping:
bowtieDir=/home-link/bcglc01/bowtie2-2.2.4

$bowtieDir/bowtie2 -p 12 -3 75 --reorder --very-sensitive -x /home-link/bcglc01/mapping/index_genome_for_bowtie2/hybrid_ColxLyrata/Ath_index ${file_R1} -S ${inputFileName}_mapped.sam

grep '^@' ${inputFileName}_mapped.sam > header_file
grep -v '^@' ${inputFileName}_mapped.sam | awk '($2==0||$2==16){print $0; next}' > filtered_file

cat header_file filtered_file > ${inputFileName}_mapped.sam

echo 'total raw reads:' > log_file_mapping_${inputFileName}.txt
LineNumber=`wc -l ${file_R1} | awk '{print $1}'`; echo $((LineNumber/4)) >> log_file_mapping_${inputFileName}.txt
echo 'total mapped reads:' >> log_file_mapping_${inputFileName}.txt
LineNumber=`wc -l filtered_file | awk '{print $1}'`; echo $((LineNumber/1)) >> log_file_mapping_${inputFileName}.txt

rm header_file filtered_file

#output as a bam file
samtools view -Sb ${inputFileName}_mapped.sam > ${inputFileName}_mapped_SE.bam
samtools sort ${inputFileName}_mapped_SE.bam -o ${inputFileName}_mapped_SE_sorted.bam
samtools index ${inputFileName}_mapped_SE_sorted.bam

rm ${inputFileName}_mapped.sam ${inputFileName}_mapped_SE.bam
