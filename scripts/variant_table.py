import pysam
import sys


def snp_sort_key(snp):
    return int(snp.contig) << 32 | snp.pos


class SNPInfo:

    def __init__(self, contig, pos, nuc, name=None, ref="C"):
        self.contig = contig
        self.pos = int(pos)
        self.nuc = nuc
        self.ref = ref
        self.name = name if name else "%s:%i:%c" % (self.contig, self.pos,
                                                    self.nuc)

    def set_allele_frequencies(self, vcfpath):
        freqs = get_allele_frequencies(vcfpath, self)
        try:
            self.quant_est_freq = freqs['quant']
            self.pass_freq = freqs['pass']
            self.raw_freq = freqs['raw']
            sys.stderr.write("Allele frequencies set for allele.\n")
        except TypeError:
            sys.stderr.write("Note: allele frequencies is null. "
                             "Allele not found.\n")
            self.quant_est_freq = None
            self.pass_freq = None
            self.raw_freq = None

    def __str__(self):
        return "%s|%s:i:c|%f:%f:f" % (self.name, self.contig, self.pos,
                                      self.nuc, self.quant_est_freq,
                                      self.pass_freq, self.raw_freq)


def make_fmt_dict(vrec):
    toks = str(vrec).strip().split("\t")
    return {k: v for k, v in zip(toks[8].split(":"), toks[9].split(":"))}


def get_allele_counts(vcfpath, snp):
    contig, pos, nuc, ref = snp.contig, snp.pos, snp.nuc, snp.ref
    assert isinstance(nuc, str) and len(nuc) == 1
    vh = pysam.VariantFile(vcfpath)
    vrec = next(vh)
    rid = vh.header.contigs.get(contig).id
    while vrec.rid != rid or vrec.pos != pos:
        vrec = next(vh)
    print(str(vrec))
    fmt = make_fmt_dict(vrec)
    try:
        index = vrec.alleles.index(nuc)
    except ValueError:
        print(vrec.alleles)
        index = -1
        sys.stderr.write("Allele %c not found. "
                         "Returning nonsense data.\n" % nuc)
    allele_passes = list(map(int, fmt["BMF_PASS"].split(",")))
    passes = list(map(int, fmt["ADP_PASS"].split(",")))
    quants = list(map(int, fmt["BMF_QUANT"].split(",")))
    quant_count = quants[index] if allele_passes[index] and index >= 0 else 0
    try:
        index_ref = vrec.alleles.index(ref)
    except ValueError:
        index_ref = -1
    quant_count_ref = quants[index_ref] if index_ref >= 0 else 0
    sys.stderr.write("alt: %i. ref: %i. sum: %i.\n" % (quant_count, quant_count_ref, sum(quants)))
    return {"alt": quant_count, "ref": quant_count_ref}


def get_allele_frequencies(vcfpath, snp):
    contig, pos, nuc = snp.contig, snp.pos, snp.nuc
    assert isinstance(nuc, str) and len(nuc) == 1
    vh = pysam.VariantFile(vcfpath)
    vrec = next(vh)
    rid = vh.header.contigs.get(contig).id
    while vrec.rid != rid or vrec.pos != pos:
        vrec = next(vh)
    print(str(vrec))
    fmt = make_fmt_dict(vrec)
    index = vrec.alleles.index(nuc) if nuc in vrec.alleles else -1
    if index == -1:
        print(vrec.alleles)
        sys.stderr.write("Allele %c not found. "
                         "Returning nonsense data.\n" % nuc)
        return {"quant": 0,
                "pass": 0.0,
                "all": 0.0,
                "quantc": 0.0,
                "quants": sum(map(int, fmt["BMF_QUANT"].split(","))),
                "passc": 0,
                "passs": sum(map(int, fmt["ADP_PASS"].split(","))),
                "allc": 0,
                "alls": sum(map(int, fmt["ADP_ALL"].split(",")))}
    allele_passes = list(map(int, fmt["BMF_PASS"].split(",")))
    passes = list(map(int, fmt["ADP_PASS"].split(",")))
    pass_count = passes[index]
    pass_sum = sum(passes)
    pass_allele_frac = pass_count * 1. / pass_sum
    all_vals = list(map(int, fmt["ADP_ALL"].split(",")))
    all_count = all_vals[index]
    all_sum = sum(all_vals)
    all_allele_frac = all_count * 1. / all_sum
    quants = list(map(int, fmt["BMF_QUANT"].split(",")))
    quant_count = quants[index] if allele_passes[index] else 0
    quant_sum = sum(quants)
    quant_est = quant_count * 1. / quant_sum
    return {"quant": quant_est,
            "pass": pass_allele_frac,
            "all": all_allele_frac,
            "quantc": quant_count,
            "quants": quant_sum,
            "passc": pass_count,
            "passs": pass_sum,
            "allc": all_count,
            "alls": all_sum}


def trim_vcfpath(vcfpath):
    return vcfpath.split(".vcf")[0]


def build_allele_count_table(outpath, vcfpaths, contig, pos, nuc):
    snp = SNPInfo(contig, pos, nuc)
    allele_count_dicts = {trim_vcfpath(vcfpath):
                          get_allele_counts(vcfpath, snp)
                          for vcfpath in vcfpaths}
    with open(outpath, "w") as f:
        f.write("#Variant: %s:%i:%c\n" % (contig, pos, nuc))
        f.write("#name\tbmf alt count\tbmf ref count\n")
        for path, acd in sorted(allele_count_dicts.items(),
                                key=lambda x: x[0] if hasattr(x, "__getitem__")
                                else 0):
            f.write("%s\t%i\t%i\t%f\n" % (path, acd['alt'], acd['ref'],
                                          acd['alt'] / (acd['alt'] +
                                                        acd['ref'])))


def build_allele_freq_table(outpath, vcfpaths, contig, pos, nuc):
    snp = SNPInfo(contig, pos, nuc)
    allele_freq_dicts = {trim_vcfpath(vcfpath):
                         get_allele_frequencies(vcfpath, snp)
                         for vcfpath in vcfpaths}
    with open(outpath, "w") as f:
        f.write("#Variant: %s:%i:%c\n" % (contig, pos, nuc))
        f.write("#name\test_freq\tpass_freq\traw_freq"
                "\test_count\test_sum"
                "\tpass_count\tpass_sum"
                "\tall_count\tall_sum\n")
        for k, v in sorted(allele_freq_dicts.items(),
                           key=lambda x: x[0] if hasattr(x, "__getitem__")
                           else 0):
            try:
                print("All count: %i. All sum: %i." % (v['allc'], v['alls']))
                f.write("%s\t%f\t%f\t%f" % (k, v['quant'], v['pass'],
                                            v['all']))
                f.write("\t%i\t%i" % (v['quantc'], v['quants']))
                f.write("\t%i\t%i" % (v['passc'], v['passs']))
                f.write("\t%i\t%i\n" % (v['allc'], v['alls']))
            except TypeError:
                f.write("%s\t0.0\t0.0\t0.0" % k)
                f.write("\t0\t0\t0\t0\t0\t0\n")


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: %s <outpath> <vcfpath> <vcfpath2> <...>\n" % sys.argv[0])
        sys.exit(1)
    contig = "7"
    pos = 55249071
    nuc = "T"
    build_allele_count_table(sys.argv[1], sys.argv[2:], contig, pos, nuc)
    sys.exit(0)
