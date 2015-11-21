import multiprocessing as mp
import subprocess
import sys

from itertools import chain

import pysam


class BMFWorkflow(object):

    def split_bam_pr(inbam, outbam):
        tmpfq = outbam.split(".bam")[0] + ".tmp.fq"
        cstr = "bam_pr -t %i -f %s %s %s" % (self.mismatch_limit, tmpfq,
                                             inbam, outbam)
        sys.stderr.write("About to call '%s'.\n", cstr)
        subprocess.check_call(cstr, shell=True)
        if(inbam not in ["-", "stdin"]):
            subprocess.check_call(["rm", inbam])
        return tmpfq


    def aln_split_rsq(self):
        cstr = "bwa mem -CTY 0 -t %i -R %s %s %s %s" % (self.threads,
                                                        self.rg_str, self.ref,
                                                        self.raw_r1,
                                                        self.raw_r2)
        if(self.ucs):
            cstr += (" | mark_unclipped -l 0 - - | bmfsort -k ucs -@ " +
                     "%i -m %s " % (self.sort_threads, self.mem_str) +
                     "{0} -T {1} -sp {1} -".format(self.mem_str, self.prefix))
        else:
            cstr += (" | bmfsort -k bmf -m %s -@ " % self.mem_str +
                     "%i -T {1} -sp {1} -".format(self.sort_threads,
                                                  self.prefix))
        sys.stderr.write("Command string: '%s'.\n", cstr)
        # Align, convert, split, and sort
        subprocess.check_call(cstr, shell=True)
        refs = pysam.FastaFile(self.ref).references
        split_bams = ["%s.split.%s.bam" % (self.prefix, ref) for ref in refs]
        tmp_bams = ["%s.split.%s.tmprsq.bam" % (self.prefix, ref) for
                    ref in refs]
        pool = mp.Pool(processes=self.threads)
        tmp_fqs = [pool.apply_async(split_bam_pr, args=(split, tmp)) for
                   split, tmp in zip(split_bams, tmp_bams)]
        if self.ucs:
            tmp_bam = "%s.tmprsq.bam"
            cstr = ("cat %s | paste -d'~' - - - - | sort " % (" ".join(tmp_fqs)) +
                    "| tr '~' '\n' | bwa mem -pCTY 0 -t %s -R " % (self.threads) +
                    "%s %s - | samtools sort -T " % (self.rg_str, self.ref) +
                     self.prefix + " -O bam -o %s -m %s " % (tmp_bam, self.mem_str))
            sys.stderr.write("Command string: '%s'.\n", cstr)
            subprocess.check_call(cstr, shell=True)
            cstr = ("cat %s %s | samtools sort -@ " % (tmp_bam, " ".join(tmp_bams)) +
                    "%i -T %s -m %s " % (self.threads, self.prefix, self.mem_str) +
                    "-o %s -" % (self.final_bam))
            sys.stderr.write("Command string: '%s'.\n", cstr)
            subprocess.check_call(cstr, shell=True)
            [subprocess.check_call("rm %s" % i for i in
                                   list(chain.from_iterable([tmp_bams,
                                                             tmp_fqs,
                                                             tmp_bam])))]
        else:
            cstr = ("cat %s | paste -d'~' - - - - | sort " % (" ".join(tmp_fqs)) +
                    "| tr '~' '\n' | bwa mem -pCTY 0 -t %s -R " % (self.threads) +
                    "%s %s - | samtools merge " % (self.rg_str, self.ref) +
                    "%s - %s" % (self.final_bam, " ".join(tmp_bams))) 
            sys.stderr.write("Command string: '%s'.\n", cstr)
            subprocess.check_call(cstr, shell=True)
            [subprocess.check_call("rm %s" % i for i in
                                   list(chain.from_iterable([tmp_bams,
                                                             tmp_fqs])))]
            
        
    def aln_rsq(self):
        cstr = "bwa mem -CTY 0 -t %i -R %s %s %s %s" % (self.threads,
                                                        self.rg_str, self.ref,
                                                        self.raw_r1,
                                                        self.raw_r2)
        if(self.ucs):
            cstr += (" | mark_unclipped -l 0 - - | bmfsort -k ucs " +
                     " -l 0 -m %s -T %s | " % (self.mem_str, self.prefix) +
                     "bam_pr -ut %i -f %s - " % (self.mismatch_limit,
                                                 self.tmp_fq) + self.tmp_bam)
        else:
            cstr += (" | bmfsort -k bmf -l 0 -m %s" % (self.mem_str) +
                     " -T %s -@ %i " % (self.prefix, self.sort_threads) +
                     "- | bam_pr -l 0 -t %i -f " % self.mismatch_limit +
                     "%s - - | samtools sort -T " % self.tmp_fq +
                     "%s.pref -O bam -m %s " % (self.prefix, self.mem_str) +
                     " -o %s -" % self.tmp_bam)
        sys.stderr.write("Command string: '%s'.\n" % cstr)
        # Align and perform rescue
        subprocess.check_call(cstr, shell=True)
        cstr = "cat %s | paste -d'~' - - - - | sort | tr '~' '\n' | "
        cstr += "bwa mem -pCTY 0 -t %i -R " % self.threads + self.rg_str
        cstr += " %s - | samtools sort -O bam -l -0 -T " % self.ref
        cstr += self.prefix + "-m %s - | samtools merge -O " % (self.mem_str)
        cstr += "bam %s %s - " % (self.final_bam, self.tmp_bam)
        sys.stderr.write("Command string: '%s'.\n" % cstr)
        # Align rescued reads, sort, merge
        subprocess.check_call(cstr, shell=True)
        subprocess.check_call("rm %s %s" % (self.tmp_bam, self.tmp_fq),
                              shell=True)

    def get_crms(self):
        cstr = "crms -dl %i -v %i -n %i -s %s " % (self.nlen, self.max_nlen,
                                                   self.prefix_len,
                                                   self.homing)
        if self.panthera:
            cstr += " -c "
        if self.gzip_output:
            cstr += " -zg %i " % self.compression
        if self.hp_threshold:
            cstr += " -t %i " % self.hp_threshold
        if self.rescaler:
            cstr += " -r %s " % self.rescaler
        cstr += (" -o %s -f " % self.tmp_prefix +
                 self.ffq_r1.split(".fastq")[0])
        cstr += " -p %i " % self.threads
        cstr += " -m %i " % self.mask

    def __init__(self):
        self.mask = None
        self.ref = None  # Path to reference
        self.ucs = None  # Flag to use unclipped start rather than pos to rescue.
        self.threads = None  # Threads to use for crms, bwa, and pigz
        self.mem_str = None  # Memory string to pass to samtools sort and bmfsort
        self.rg_str = None  # RG string to pass to bwa
        self.opts = None  # bwa alignment options
        self.mismatch_limit = None  # Mismatch limit for rescue step.
        self.prefix = None  # temporary prefix for temporary files
        self.split = None  # Whether or not to split the bam and apply parallel rescue
        self.raw_r1 = None
        self.raw_r2 = None
        self.index_read = None
        self.ffq_r1 = None
        self.ffq_r2 = None
        self.chem = None
        self.exp = None
        self.sort_threads = None
        self.tmp_fq = None  # Used for temporary fastq file for single rescue
        self.tmp_bam = None
        self.final_bam = None
        self.homing = None


    def fq_dmp(self):
        if self.chem == 0:  # annealed
            subprocess.check_call(self.get_fqmsi(), shell=True)
        elif self.chem == 1: # loeb
            subprocess.check_call(self.get_crms(), shell=True)
        elif self.chem == 2:  # secondary index
            subprocess.check_call(self.get_fqms(), shell=True)
        else:
            raise ValueError("Chemistry not supported. Abort!")

    def align_rsq(self):
        if(self.split):
            self.split_aln_rsq()
        else:
            self.aln_rsq()



'''
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
prefix=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)
memstr="12G"
sort_threads="4"
mismatch_limit="2"


# Run-level parameters
r1=$1
r2=$2
tmp_fq=$3
final_bam=$4

# Check to see if temporary fastq path set.
# If not, set to a trimmed version of R1
[ -z "$3" ] && tmp_fq=${r1:0:`expr index $r1 fastq` - 1}.tmp.fq
[ -z "$4" ] && final_bam=${r1:0:`expr index $r1 fastq` - 1}.rsq.sort.bam

## BIG PIPE
if [ $MT_RESCUE -eq 0 ]
then
    # Align, sort, rescue.
    tmp_bam=${r1:0:`expr index $r1 fastq` - 1}.tmp.sort.bam
    echo bwa mem $opts $r1 $r2
    if [ $USE_UCS -eq 0]
    then
        bwa mem $opts $r1 $r2 | \
        bmfsort -T $prefix -m $memstr -@ $sort_threads -k bmf -l 0 | \
        bam_pr -f $tmp_fq -t $mismatch_limit - $tmp_bam
    else
        bwa mem $opts $r1 $r2 | mark_unclipped -l 0 - - \
        bmfsort -T $prefix -m $memstr -@ $sort_threads -k ucs -l 0 | \
        bam_pr -uf $tmp_fq -t $mismatch_limit - $tmp_bam
    fi

    # Cat temporary fastq, sort by read name, pipe to bwa in interleaved mode
    # and merge with original bam
    if [ $USE_UCS -eq 0]
    then
        cat $tmp_fq | paste -d'~' - - - - | sort | tr '~' '\n' | \
        bwa mem -p $opts - | samtools merge -@ $sort_threads $final_bam $tmp_bam -
        # Clean up
        rm $tmp_fq $tmp_bam
    else
        # Sort, align, sort, and convert tmp_fq
        tmprsqbam=${tmp_fq}.tmprsq.bam
        cat $tmp_fq | paste -d'~' - - - - | sort | tr '~' '\n' | \
        bwa mem -p $opts - | \
        samtools sort -O bam -T $prefix -m $memstr -o $tmprsqbam

        # Sort tmp_bam and merge in with tmprsqbam
        samtools sort -O bam -T $prefix -m $memstr $tmprsqbam | \
        samtools merge -O bam -l 0 -@ $sort_threads $final_bam $tmp_bam -

        # Clean up
        rm $tmprsqbam $tmp_bam $tmp_fq
    fi        
    # Index
    samtools index $final_bam

### Moderate pipe using pos
else
    tmp_bam=${r1:0:`expr index $r1 fastq` - 1}.tmp.sort.bam
    split_prefix=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)

    # Align
    if [ $USE_UCS -eq 0 ]
    then
        bwa mem $opts $r1 $r2 | \
        bmfsort -T $prefix -m $memstr -@ $sort_threads -k bmf -sp $split_prefix
    else
        bwa mem $opts $r1 $r2 | mark_unclipped -l 0 - - | \
        bmfsort -T $prefix -m $memstr -@ $sort_threads -k ucs -sp $split_prefix
    fi

    # Apply bam rescue in parallel

    
    if [ $USE_UCS -eq 0 ]
    then
        for tmp in $(ls $split_prefix*.bam)
        do
            sem -j $threads bam_pr -f ${tmp%.bam}.tmp.fq -t $mismatch_limit $tmp ${tmp%.bam}.rsq.bam
        done
        sem --wait
    else
        for tmp in $(ls $split_prefix*.bam)
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
        rm $(ls ${split_prefix}*.bam ${split_prefix}*.tmp.fq $tmp_bam)
    else
        tmprsqbam=${split_prefix}.dnd.bam
        # Sort the fastq, pipe to bwa, pipe to big merge.
        cat $(ls ${split_prefix}*.tmp.fq) | \
        # Put all tmp fastqs into a stream and sort by read name
        paste -d'~' - - - - | sort | tr '~' '\n' | \
        # Align
        bwa mem -p $opts - | samtools sort -@ $sort_threads -m $memstr -T $prefix -o $tmprsqbam

        # Sort and merge bams
        cat $tmprsq $(ls ${split_prefix}*.rsq.bam) | \
        samtools sort -m $memstr -T $prefix -O bam -@ sort_threads -o $final_bam

        # Clean up
        rm $(ls ${split_prefix}*.bam ${split_prefix}*.tmp.fq $tmp_bam $tmprsqbam)
    fi
    # Index bam
    samtools index $final_bam

fi

echo "Successfully produced $final_bam"
return 0
    

    

'''