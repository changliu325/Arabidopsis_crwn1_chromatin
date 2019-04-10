#!/bin/bash

# !!!!Replace it with your own directory!!!
SICER=/Users/chang/Documents/ZMBP/tools_and_resources/SICER_V1.1/SICER
Data_source=/Users/chang/Documents/ZMBP/DNA_methylation_and_NP/ChIP-seq

#["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"]   ["gap size (bp)"] ["FDR"]
sh $SICER/SICER.sh $Data_source CRWN1_2HA_11_3_antiHA_rep1_R1.fastq_mapped_PE_sorted.bam.bed CRWN1_2HA_11_3_input_rep1_R1.fastq_mapped_PE_sorted.bam.bed . tair10 5 500 300 0.88 1500 .01

sh $SICER/SICER.sh $Data_source CRWN1_2HA_4_4_antiHA_R1.fastq_mapped_PE_sorted.bam.bed CRWN1_2HA_4_4_input_R1.fastq_mapped_PE_sorted.bam.bed . tair10 5 500 300 0.88 1500 .01
