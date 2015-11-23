import multiprocessing as mp
import subprocess
import sys
import uuid

from enum import Enum
from itertools import chain

import pysam

class chemistry(Enum):
    loeb = 0
    annealed = 1
    secondary = 2

def split_fq(instr):
    return instr.split(".fastq")[0].split(".fq")[0]


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
        cstr = "crms -dl %i -v %i -s %s " % (self.blen, self.max_blen,
                                             self.homing)
        if self.panthera:
            cstr += " -c "
        if self.gzip_output:
            cstr += " -zg %i " % self.compression
        if self.hp_threshold:
            cstr += " -t %i " % self.hp_threshold
        if self.rescaler:
            cstr += " -r %s " % self.rescaler
        if self.nucsplit_n:
            cstr += " -n %i" % self.nucsplit_n
        cstr += (" -o %s -f " % self.tmp_prefix +
                 self.ffq_r1.split(".fastq")[0])
        cstr += " -p %i " % self.threads
        cstr += " -m %i " % self.mask

    def get_params(self, path):
        for line in open(path, "r").readlines():
            key, value = line.strip().split("=")  # Unsafe
            exec("self.%s = '%s'" % (key, value))

    def check_params(self):
        self.ffq_r1 = split_fq(self.raw_r1) + ".r1.dmp.fq"
        self.ffq_r2 = split_fq(self.raw_r2) + ".r2.dmp.fq"
        self.final_bam = split_fq(self.raw_r1) + ".final.bam"

        self.mismatch_limit = int(self.mismatch_limit) if self.mismatch_limit else 2

        self.gzip_output = (self.gzip_output.lower() == 'true' or
                            int(self.gzip_output)) if self.gzip_output else True

        self.compression = int(self.compression) if self.compression else 6

        self.mask = int(self.mask) if self.mask else 1

        self.hp_threshold = int(self.hp_threshold) if self.hp_threshold else 12

        if self.blen:
            self.blen = int(self.blen)
        elif self.chem in [chemistry.loeb, chemistry.annealed]:
            raise ValueError("barcode length required for inline chemistry.")

        self.max_blen = int(self.max_blen) if self.max_blen else self.blen
        if not self.mem_str:
            self.mem_str = "3G"
            sys.stderr.write("mem_str unset. Setting to 3G.\n")

        if self.panthera:
            self.panthera = 1
            sys.stderr.write("mem_str unset. Setting to 3G.\n")
        else:
            self.panthera = (self.panthera.lower() == "true" or
                             int(settings.panthera))

        if self.prefix is None:
            self.prefix = uuid.uuid4().get_hex().upper()[:16]
            sys.stderr.write("Prefix unset. Random: %s.\n" % self.prefix)

        if self.nucsplit_n:
            self.nucsplit_n = int(self.nucsplit_n)

        if self.ref is None:
            raise ValueError("Required reference path missing!")

        if self.rg_str is None:
            self.rg_str = "@RG\tID:omgz\tSM:wtf\tPL:ILMN\tDS:%s" % self.prefix
        sys.stderr.write("RG String: %s.\n" % self.rg_str)

        self.sort_threads = int(self.sort_threads) if self.sort_threads else 2

        self.split = (self.split.lower() == "true" or
                      int(self.split)) if self.split else False

        if self.threads:
            self.threads = int(self.threads)
        else:
            self.threads = 4
            sys.stderr.write("threads unset. Defaulting to 4.\n")

        if not self.tmp_bam:
            self.tmp_bam = self.prefix + ".tmp.bam"

        if not self.tmp_fq:
            self.tmp_fq = self.prefix + ".tmp.fq"

        self.ucs = (self.ucs.lower() == "true" or
                    int(self.ucs)) if(self.ucs) else False
        '''
        if self.ucs:
            self.ucs = (self.ucs.lower() == "true" or int(self.ucs))
        else:
            self.ucs = False
        '''

    def __init__(self, r1, r2, rescaler, index, adapter_design):
        if adapter_design == "l":
            self.chem = chemistry.loeb
        elif adapter_design == "a":
            self.chem = chemistry.annealed
        elif adapter_design == "s":
            self.chem = chemistry.secondary
        else:
            raise ValueError("Not valid chemistry. Choose 'l', 'a', or 's'.")
        if self.chem == chemistry.secondary and index in None:
            raise ValueError("For secondary index chemistry, an index fastq must be provided.")
        self.compression = None
        self.ffq_r1 = None
        self.ffq_r2 = None
        self.final_bam = None
        self.gzip_output = None
        self.homing = None
        self.hp_threshold = None
        self.index_read = index_read
        self.mask = None
        self.max_blen = None
        self.mem_str = None  # Memory string to pass to samtools sort and bmfsort
        self.mismatch_limit = None  # Mismatch limit for rescue step.
        self.blen = None
        self.panthera = 1
        self.prefix = None  # temporary prefix for temporary files
        self.nucsplit_n = None
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
        if self.chem == chemistry.loeb:
            subprocess.check_call(self.get_crms(), shell=True)
        elif self.chem == chemistry.annealed:
            subprocess.check_call(self.get_fqmsi(), shell=True)
        elif self.chem == chemistry.secondary:
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
    parser.add_argument("--index", "-i",
                        help="Path to index read. Leave unset for inline barcodes.")
    parser.add_argument("--conf", "-c", help="Path to config file.")
    parser.add_argument("--adapter-design", "-a", help="Adapter design. "
                        "'l' for loeb,'a' for annealed, 's' for secondary.",
                        required=True)
    return parser.parse_args()

def main():
    args = get_args()
    wf = BMFWorkflow(args.r1, args.r2, args.rescaler, args.index, args.adapter_design)
    wf.get_params(args.conf)  # Parse parameters from config file
    wf.check_params()  # Make sure that all required parameters are set
    wf.fq_dmp()  # Perform molecular demultiplexing
    wf.aln_rsq()
    return 0

if __name__ == "__main__":
    sys.exit(main())