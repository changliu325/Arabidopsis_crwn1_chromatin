input=$1
output=$2

#format fastq file, so every read will have 101bp

awk '(NR % 4==2){seq=$0; Ns="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; all=(seq""Ns); print substr(all,1,101); next}(NR % 4==0){seq=$0; Ns="#####################################################################################################"; all=(seq""Ns); print substr(all,1,101); next}{print $0}' $input > $output
