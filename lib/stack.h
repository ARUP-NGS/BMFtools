#ifndef UNIQUE_OBS_H
#define UNIQUE_OBS_H

#include "dlib/bam_util.h"
#include "dlib/vcf_util.h"
#include <tgmath.h>

#include <unordered_map>

#define DEFAULT_MAX_DEPTH (1 << 18)

namespace BMF {

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
        "##INFO=<ID=SOMATIC_CALL,Number=R,Type=Integer,Description=\"Boolean value for a somatic call for each allele.\">",
        "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">"
};

    struct stack_conf_t {
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
    uint32_t agreed;
    uint32_t size;
    public:
        int is_pass() {
            return pass;
        }
        uint32_t get_size() {
            return size;
        }
        uint32_t get_agreed() {
            return agreed;
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
            return is_duplex1 + (is_duplex2 >= 0 ? is_duplex2: 0);
        }
        UniqueObservation() {
            memset(this, 0, sizeof(*this));
        }
        UniqueObservation(const bam_pileup1_t& plp):
            qname(bam_get_qname(plp.b)),
            cycle1(dlib::arr_qpos(&plp)),
            cycle2(-1),
            discordant(-1),
            quality(((uint32_t *)dlib::array_tag(plp.b, "PV"))[cycle1]),
            mq1(plp.b->core.qual),
            mq2((uint32_t)-1),
            rv((uint32_t)bam_itag(plp.b, "RV")),
            is_duplex1(bam_itag(plp.b, "DR")),
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
            agreed(((uint32_t *)dlib::array_tag(plp.b, "FA"))[cycle1]),
            size(bam_itag(plp.b, "FM"))
        {
        }
        void add_obs(const bam_pileup1_t& plp);
    };

    class stack_aux_t {
    public:
        stack_conf_t conf;
        dlib::BamHandle tumor;
        dlib::BamHandle normal;
        dlib::VcfHandle vcf;
        faidx_t *fai;
        khash_t(bed) *bed;
        int last_tid;
        char *ref_seq;
        void pair_region_itr(int bamtid, int start, int stop, int &tn_plp, int &tpos, int &ttid, int &nn_plp, int &npos, int &ntid) {
            bam_plp_reset(tumor.plp);
            bam_plp_reset(normal.plp);
            if(tumor.iter) hts_itr_destroy(tumor.iter);
            if(normal.iter) hts_itr_destroy(normal.iter);
            tumor.iter = bam_itr_queryi(tumor.idx, bamtid, start, stop);
            normal.iter = bam_itr_queryi(normal.idx, bamtid, start, stop);
            while((tumor.pileups = bam_plp_auto(tumor.plp, &ttid, &tpos, &tn_plp)) != 0) {
                if(tpos < start && ttid == bamtid) continue;
                if(tpos >= start) {
                    //LOG_DEBUG("Breaking? %i >= start (%i)\n", tpos, start);
                    break;
                }
                LOG_EXIT("Wrong tid (ttid: %i, bamtid %i)? wrong pos? tpos, stop %i, %i", ttid, bamtid, tpos, stop);
            }
            while((normal.pileups = bam_plp_auto(normal.plp, &ntid, &npos, &nn_plp)) != 0) {
                if(npos < start && ntid == bamtid) continue;
                if(npos >= start) {
                    //LOG_DEBUG("Breaking? %i >= start (%i)\n", npos, start);
                    break;
                }
            }
            assert(npos == tpos && ntid == ttid);
            //LOG_DEBUG("Advanced to %i:%i/%i:%i.\n", ttid, bamtid, npos, start);
            //LOG_DEBUG("N pileups: %i, %i.\n", tn_plp, nn_plp);
        }
        int next_paired_pileup(int *ttid, int *tpos, int *tn_plp, int *ntid, int *npos, int *nn_plp, int stop) {
            if((tumor.pileups = bam_plp_auto(tumor.plp, ttid, tpos, tn_plp)) != 0 &&
                    (normal.pileups = bam_plp_auto(normal.plp, ntid, npos, nn_plp)) != 0) {
                return *npos < stop;
            }
            return 0;
        }
        stack_aux_t(char *tumor_path, char *normal_path, char *vcf_path, bcf_hdr_t *vh, stack_conf_t conf):
            conf(conf),
            tumor(tumor_path),
            normal(normal_path),
            vcf(vcf_path, vh, conf.output_bcf ? "wb": "w"),
            fai(nullptr),
            bed(nullptr),
            last_tid(-1),
            ref_seq(nullptr)
        {
            dlib::bcf_add_bam_contigs(vcf.vh, tumor.header);
            if(!conf.max_depth) conf.max_depth = DEFAULT_MAX_DEPTH;
        }
        char get_ref_base(int tid, int pos) {
            //LOG_DEBUG("fai ptr %p.\n", (void *)fai);
            int len;
            if(tid != last_tid) {
                //LOG_DEBUG("Loading ref_seq\n");
                if(ref_seq) free(ref_seq);
                ref_seq = fai_fetch(fai, tumor.header->target_name[tid], &len);
                last_tid = tid;
            }
            //LOG_DEBUG("Returning ref base\n");
            return ref_seq[pos];
        }
        ~stack_aux_t() {
            if(bed) dlib::bed_destroy_hash((void *)bed);
            if(ref_seq) free(ref_seq);
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
                        normal(nobs, _tid, _pos)
        {
        }
    };
}

#endif
