#!/bin/bash
# ml load dtu_polisher/1.0 kma/1.4.4 minimap2/2.17 miniasm/0.3r179 samtools/1.9

# Define the options and their default values
taxadb=""
resistancedb=""
input_file=""
output_file=""
num_cores=$(nproc)
threads=$((num_cores/2))

# Use getopts to parse the options
while getopts ":t:r:i:o:|threads:" opt; do
  case $opt in
    t) input_file=$OPTARG
       ;;
    r) output_file=$OPTARG
       ;;
    i) input_file=$OPTARG
       ;;
    o) output_file=$OPTARG
       ;;
    threads) threads=$OPTARG
       ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :) echo "Option -$OPTARG requires an argument." >&2
       exit 1
       ;;
  esac
done

# Print usage documentation if no options are passed
if [[ -z "$taxadb" ]] || [[ -z "$resistancedb" ]] || [[ -z "$input_file" ]] || [[ -z "$output_file" ]]; then
  echo "Usage: myscript.sh -t taxa_dabase -r resistance_database -i input_file -o output_file"
  exit 1
fi

# run KMA: taxa
kma -i $input_file -o ${output_file}_taxa -t_db $taxadb -tmp $(pwd) \
-mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nc -nf -1t1 -ca -verbose 2 -t $threads \
-sam > ${output_file}_taxa.sam

# number of samtools threads (-@) needs to be the same as used for KMA. Otherwise risk of segmentation fault
echo "converting sam to bam"
samtools view -u -h -@ 45 -o ${output_file}_taxa.bam ${output_file}_taxa.sam 
rm ${output_file}_taxa.sam 

echo "Sorting bam"
samtools sort -@ 45 -o ${output_file}_taxa_sort.bam ${output_file}_taxa.bam
rm ${output_file}_taxa.bam
mv ${output_file}_taxa_sort.bam ${output_file}_taxa.bam

# run KMA: resistance
kma -i $input_file -o ${output_file}_ARG -t_db $resistancedb -tmp $(pwd) \
-mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nc -nf -1t1 -ca -verbose 2 -t $threads \
-sam > ${output_file}_ARG.sam

# number of samtools threads (-@) needs to be the same as used for KMA. Otherwise risk of segmentation fault
echo "converting sam to bam"
samtools view -u -h -@ 45 -o ${output_file}_ARG.bam ${output_file}_ARG.sam 
rm ${output_file}_ARG.sam 

echo "Sorting bam"
samtools sort -@ 45 -o ${output_file}_ARG_sort.bam ${output_file}_ARG.bam
rm ${output_file}_ARG.bam
mv ${output_file}_ARG_sort.bam ${output_file}_ARG.bam