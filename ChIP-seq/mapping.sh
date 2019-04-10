file_R1=$1
file_R2=$2
#extract the file name (in case directory is included in input)
inputFileName=${file_R1##*/}
#input files, the raw fastq files of read1 

#mapping:

bowtieDir=/home-link/bcglc01/bowtie2-2.2.4

$bowtieDir/bowtie2 -p 12 --local --reorder --very-sensitive -x /home-link/bcglc01/mapping/index_genome_for_bowtie2/hybrid_ColxLyrata/Ath_index -1 $file_R1 -2 $file_R2 -S ${inputFileName}_mapped.sam

#filter sam files: 
# a) both ends mapped, insert size shall not be larger than 1.5kb

grep '^@' ${inputFileName}_mapped.sam > header_file
grep -v '^@' ${inputFileName}_mapped.sam | awk '($9==0||$9>1500||$9<-1500){next}{print $0}' > filtered_file

cat header_file filtered_file > ${inputFileName}_mapped.sam

#generate a temporal file to check sequence coverage, each line tells you the mapping chromosome and position of read1:

awk '(NR % 2==1){print $3"\t"$4;next}{next}' filtered_file > ${inputFileName}_coverage_estimate.txt

rm header_file filtered_file

#output as a bam file
samtools view -Sb ${inputFileName}_mapped.sam > ${inputFileName}_mapped_PE.bam
samtools sort ${inputFileName}_mapped_PE.bam -o ${inputFileName}_mapped_PE_sorted.bam
samtools index ${inputFileName}_mapped_PE_sorted.bam

rm ${inputFileName}_mapped.sam ${inputFileName}_mapped_PE.bam
