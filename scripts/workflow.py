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

    def aln_pipe_rsq(self):
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

    def get_params(self, path):
        for line in open(path, "r").readlines():
            key, value = line.strip().split("=")  # Unsafe
            exec("self.%s = '%s'" % (key, value))

    def __init__(self, r1, r2, rescaler):
        self.chem = None
        self.compression = None
        self.ffq_r1 = None
        self.ffq_r2 = None
        self.final_bam = None
        self.gzip_output = None
        self.homing = None
        self.hp_threshold = None
        self.index_read = None
        self.mask = None
        self.max_nlen = None
        self.mem_str = None  # Memory string to pass to samtools sort and bmfsort
        self.mismatch_limit = None  # Mismatch limit for rescue step.
        self.nlen = None
        self.panthera = None
        self.prefix = None  # temporary prefix for temporary files
        self.prefix_len = None
        self.raw_r1 = r1
        self.raw_r2 = r2
        self.rescaler = rescaler
        self.ref = None  # Path to reference
        self.rg_str = None  # RG string to pass to bwa
        self.sort_threads = None
        self.split = None  # Whether or not to split the bam and apply parallel rescue
        self.threads = None  # Threads to use for crms, bwa, and pigz
        self.tmp_bam = None
        self.tmp_fq = None  # Used for temporary fastq file for single rescue
        self.ucs = None  # Flag to use unclipped start rather than pos to rescue.


    def fq_dmp(self):
        if self.chem == 0:  # annealed
            subprocess.check_call(self.get_fqmsi(), shell=True)
        elif self.chem == 1: # loeb
            subprocess.check_call(self.get_crms(), shell=True)
        elif self.chem == 2:  # secondary index
            subprocess.check_call(self.get_fqms(), shell=True)
        else:
            raise ValueError("Chemistry not supported. Abort!")

    def aln_rsq(self):
        if(self.split):
            self.aln_split_rsq()
        else:
            self.aln_pipe_rsq()

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("r1", help="Path to raw r1 fastq.")
    parser.add_argument("r2", help="Path to raw r2 fastq.")
    parser.add_argument("--rescaler", "-r",
                        help="Path to rescaler array in human-readable format.")
    parser.add_argument("--conf", "-c", help="Path to config file.")
    return parser.parse_args()

def main():
    args = get_args()
    wf = BMFWorkflow(args.r1, args.r2, args.rescaler)
    wf.get_params(args.conf)
    wf.fq_dmp()
    wf.aln_rsq()
    return 0

if __name__ == "__main__":
    sys.exit(sys.exit(0))