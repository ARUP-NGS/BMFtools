#!/bin/bash

# Exit upon failure
set -e
# Print all statements
set -x

# Parameters
USE_UCS=1
MT_RESCUE=1
ref="/mounts/genome/human_g1k_v37.fasta"
threads="8"
rg_str="@RG\tID:omgz\tSM:wtf\tPL:ILMN"
opts="-CYT 0 -t $threads -v 1 -R $rg_str $ref "
tmp_prefix=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)
memstr="12G"
sort_threads="4"
mismatch_limit="2"


# Run-level parameters
r1=$1
r2=$2
tmpfq=$3
final_bam=$4

# Check to see if temporary fastq path set.
# If not, set to a trimmed version of R1
[ -z "$3" ] && tmpfq=${r1:0:`expr index $r1 fastq` - 1}.tmp.fq
[ -z "$4" ] && final_bam=${r1:0:`expr index $r1 fastq` - 1}.rsq.sort.bam

## BIG PIPE
if [ $MT_RESCUE -eq 0 ]
then
	# Align, sort, rescue.
	tmpbam=${r1:0:`expr index $r1 fastq` - 1}.tmp.sort.bam
	echo bwa mem $opts $r1 $r2
	if [ $USE_UCS -eq 0]
	then
		bwa mem $opts $r1 $r2 | \
		bmfsort -T $tmp_prefix -m $memstr -@ $sort_threads -k bmf -l 0 | \
		bam_pr -f $tmpfq -t $mismatch_limit - $tmpbam
	else
		bwa mem $opts $r1 $r2 | mark_unclipped -l 0 - - \
		bmfsort -T $tmp_prefix -m $memstr -@ $sort_threads -k ucs -l 0 | \
		bam_pr -uf $tmpfq -t $mismatch_limit - $tmpbam
	fi

	# Cat temporary fastq, sort by read name, pipe to bwa in interleaved mode
	# and merge with original bam
	if [ $USE_UCS -eq 0]
	then
		cat $tmpfq | paste -d'~' - - - - | sort | tr '~' '\n' | \
		bwa mem -p $opts - | samtools merge -@ $sort_threads $final_bam $tmpbam -
		# Clean up
		rm $tmpfq $tmpbam
	else
		# Sort, align, sort, and convert tmpfq
		tmprsqbam=${tmpfq}.tmprsq.bam
		cat $tmpfq | paste -d'~' - - - - | sort | tr '~' '\n' | \
		bwa mem -p $opts - | \
		samtools sort -O bam -T $tmp_prefix -m $memstr -o $tmprsqbam

		# Sort tmpbam and merge in with tmprsqbam
		samtools sort -O bam -T $tmp_prefix -m $memstr $tmprsqbam | \
		samtools merge -O bam -@ $sort_threads $final_bam $tmpbam -

		# Clean up
		rm $tmprsqbam $tmpbam $tmpfq
	fi		
	# Index
	samtools index $final_bam

### Moderate pipe using pos
else
	tmpbam=${r1:0:`expr index $r1 fastq` - 1}.tmp.sort.bam
	split_prefix=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)

	# Align
	if [ $USE_UCS -eq 0 ]
	then
		bwa mem $opts $r1 $r2 | \
		bmfsort -T $tmp_prefix -m $memstr -@ $sort_threads -k bmf -sp $split_prefix
	else
		bwa mem $opts $r1 $r2 | mark_unclipped -l 0 - - | \
		bmfsort -T $tmp_prefix -m $memstr -@ $sort_threads -k ucs -sp $split_prefix
	fi

	# Apply bam rescue in parallel

	
	if [ $USE_UCS -eq 0 ]
	then
		for tmp in $(ls ${split_prefix}*.bam)
		do
			sem -j $threads bam_pr -f ${tmp%.bam}.tmp.fq -t $mismatch_limit $tmp ${tmp%.bam}.rsq.bam
		done
		sem --wait
	else
		for tmp in $(ls ${split_prefix}*.bam)
		do
			sem -j $threads bam_pr -uf ${tmp%.bam}.tmp.fq -t $mismatch_limit $tmp ${tmp%.bam}.rsq.bam
		done
		sem --wait
	fi

	if [ $USE_UCS -eq 0 ]
	then
		# Sort the fastq, pipe to bwa, pipe to big merge.
		cat $(ls ${split_prefix}*.tmp.fq) | \
		# Put all tmp fastqs into a stream and sort by read name
		paste -d'~' - - - - | sort | tr '~' '\n' | \
		# Align
		bwa mem -p $opts - | samtools merge -@ $sort_threads $final_bam $(ls ${split_prefix}*.rsq.bam)
		# rm $(ls ${split_prefix}*.bam ${split_prefix}*.tmp.fq $tmpbam)
	else
		tmprsqbam=${split_prefix}.dnd.bam
		# Sort the fastq, pipe to bwa, pipe to big merge.
		cat $(ls ${split_prefix}*.tmp.fq) | \
		# Put all tmp fastqs into a stream and sort by read name
		paste -d'~' - - - - | sort | tr '~' '\n' | \
		# Align
		bwa mem -p $opts - | samtools sort -@ $sort_threads -m $memstr -T $tmp_prefix -o $tmprsqbam

		# Sort and merge bams
		cat $tmprsq $(ls ${split_prefix}*.rsq.bam) | \
		samtools sort -m $memstr -T $tmp_prefix -O bam -@ sort_threads -o $final_bam

		# Clean up
		# rm $(ls ${split_prefix}*.bam ${split_prefix}*.tmp.fq $tmpbam $tmprsqbam)
	fi
	# Index bam
	samtools index $final_bam

fi

echo "Successfully produced $final_bam"
return 0
	

	
