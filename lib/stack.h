#ifndef UNIQUE_OBS_H
#define UNIQUE_OBS_H

#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include "include/igamc_cephes.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "dlib/bed_util.h"
#include "dlib/bam_util.h"
#include "dlib/misc_util.h"
#include "dlib/vcf_util.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#define DEFAULT_MAX_DEPTH (1 << 18)

static const char *stack_vcf_lines[] = {
        "##FORMAT=<ID=BMF_PASS,Number=R,Type=Integer,Description=\"1 if variant passes, 0 otherwise.\">",
        "##FORMAT=<ID=ADP,Number=R,Type=Integer,Description=\"Number of unique observations for each allele.\">",
        "##FORMAT=<ID=ADPO,Number=R,Type=Integer,Description=\"Number of unique observations of overlapped read pairs for each allele.\">",
        "##FORMAT=<ID=ADPD,Number=R,Type=Integer,Description=\"Number of duplex observations for each allele. If both reads in an overlapping pair are duplex, this counts each separately.\">",
        "##FORMAT=<ID=ADPR,Number=R,Type=Integer,Description=\"Total number of original reversed reads supporting allele.\">",
        "##FORMAT=<ID=RVF,Number=R,Type=Float,Description=\"Fraction of reads supporting allele which were reversed.\">",
        "##FORMAT=<ID=QSS,Number=R,Type=Integer,Description=\"Q Score Sum for each allele for each sample.\">",
        "##FORMAT=<ID=AMBIG,Number=1,Type=Integer,Description=\"Number of ambiguous (N) base calls at position.\">",
        "##INFO=<ID=SOMATIC_PV,Number=R,Type=Float,Description=\"P value for a somatic call for each allele.\">",
        "##INFO=<ID=SOMATIC_CALL,Number=R,Type=Integer,Description=\"P value for a somatic call for each allele.\">"
};



namespace BMF {

    struct stack_conf {
        float minFR; // Minimum fraction of family members agreed on base
        float minAF; // Minimum aligned fraction
        int max_depth;
        uint32_t minFM;
        uint32_t minFA;
        uint32_t minPV;
        uint32_t minMQ;
        int minCount;
        int minDuplex;
        int minOverlap;
        int skip_improper;
        uint32_t skip_flag; // Skip reads with any bits set to true
        int output_bcf;
    };

    class SampleVCFPos;

    class UniqueObservation {
    friend SampleVCFPos;
    std::string qname;
    int cycle1;
    int cycle2; // Masked, from other read, if it was found.
    int discordant;
    uint32_t quality;
    uint32_t mq1;
    uint32_t mq2;
    uint32_t rv;
    int is_duplex1;
    int is_duplex2;
    int is_reverse1;
    int is_reverse2;
    int is_overlap;
    int pass;
    double pvalue;
    int flag; // May turn into a flag
    char base1;
    char base2; // Masked, from other read
    char base_call;
    public:
        uint32_t agreed;
        uint32_t size;
        int is_pass() {
            return pass;
        }
        int get_reverse() {
            return is_reverse1 + (is_reverse2 >= 0) ? is_reverse2: 0;
        }
        void set_pass(int _pass) {
            pass = _pass;
        }
        uint32_t get_meanMQ() {
            return mq2 == (uint32_t)-1 ? mq1: ((mq2 + mq1 + 0.5) / 2);
        }
        double get_FA() {
            return (double)agreed / size;
        }
        int get_overlap() {return is_overlap;}
        uint32_t get_quality() {return quality;}
        int get_duplex() {
            return is_duplex1 + (is_duplex2 != -1 ? is_duplex2: 0);
        }
        UniqueObservation() {
            memset(this, 0, sizeof(*this));
        }
        inline UniqueObservation(const bam_pileup1_t& plp):
            qname(bam_get_qname(plp.b)),
            cycle1(arr_qpos(&plp)),
            cycle2(-1),
            discordant(-1),
            quality(((uint32_t *)array_tag(plp.b, "PV"))[cycle1]),
            mq1(plp.b->core.qual),
            mq2((uint32_t)-1),
            rv((uint32_t)bam_aux2i(bam_aux_get(plp.b, "RV"))),
            is_duplex1(bam_aux2i(bam_aux_get(plp.b, "DR"))),
            is_duplex2(-1),
            is_reverse1((plp.b->core.flag & BAM_FREVERSE) != 0),
            is_reverse2(-1),
            is_overlap(0),
            pass(1),
            pvalue(std::pow(10, quality - 0.1)),
            flag(plp.b->core.flag),
            base1(seq_nt16_str[bam_seqi(bam_get_seq(plp.b), plp.qpos)]),
            base2('\0'),
            base_call(base1),
            agreed(((uint32_t *)array_tag(plp.b, "FA"))[cycle1]),
            size(bam_aux2i(bam_aux_get(plp.b, "FM")))
        {
        }
        void add_obs(const bam_pileup1_t& plp);
    };

    class stack_aux_t {
    public:
        dlib::BamHandle *tumor;
        dlib::BamHandle *normal;
        dlib::VcfHandle *vcf;
        faidx_t *fai;
        khash_t(bed) *bed;
        int last_tid;
        char *ref_seq;
        stack_conf conf;
        stack_aux_t(char *tumor_path, char *normal_path, char *vcf_path, bcf_hdr_t *vh, BMF::stack_conf _conf) {
            memset(this, 0, sizeof(*this));
            last_tid = -1;
            tumor = new dlib::BamHandle(tumor_path);
            normal = new dlib::BamHandle(normal_path);
            vcf = new dlib::VcfHandle(vcf_path, vh, _conf.output_bcf ? "wb": "w");
            bcf_add_bam_contigs(vcf->vh, tumor->header);
            conf = _conf;
            if(!conf.max_depth) conf.max_depth = DEFAULT_MAX_DEPTH;
        }
        char get_ref_base(int tid, int pos) {
            int len;
            if(tid != last_tid) {
                if(ref_seq) free(ref_seq);
                ref_seq = fai_fetch(fai, tumor->header->target_name[tid], &len);
                last_tid = tid;
            }
            return ref_seq[pos];
        }
        ~stack_aux_t() {
            if(bed) bed_destroy_hash((void *)bed);
            if(ref_seq) free(ref_seq);
            if(tumor) delete tumor;
            if(normal) delete normal;
            if(vcf) delete vcf;
        }
    };

    class PairVCFPos;

    class SampleVCFPos {
        friend PairVCFPos;
        std::unordered_map<char, std::vector<UniqueObservation *>> templates;
        size_t size;
        int32_t pos;
        int32_t tid;
    public:
        void to_bcf(bcf1_t *vrec, bcf_hdr_t *hdr, char refbase);
        SampleVCFPos(std::unordered_map<std::string, UniqueObservation>& obs, int32_t _tid, int32_t _pos);
    };

    class PairVCFPos {
        SampleVCFPos tumor;
        SampleVCFPos normal;
    public:
        void to_bcf(bcf1_t *vrec, stack_aux_t *aux, int ttid, int tpos);
        PairVCFPos(std::unordered_map<std::string, UniqueObservation>& tobs,
                std::unordered_map<std::string, UniqueObservation>& nobs,
                    int32_t _tid, int32_t _pos):
                        tumor(tobs, _tid, _pos),
                        normal(nobs, _tid, _pos) {
        }
    };
}

#endif
