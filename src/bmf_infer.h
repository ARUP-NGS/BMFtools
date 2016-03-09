#ifndef BMF_INFER_H
#define BMF_INFER_H
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <tgmath.h>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include "htslib/sam.h"
#include "include/sam_opts.h"
#include "include/bam.h" // for bam_get_library
#include "include/igamc_cephes.h" // for igamc
#include "lib/kingfisher.h"
#include "dlib/cstr_util.h"
#include "dlib/sort_util.h"
#include "dlib/bam_util.h"
#include "bmf_rsq.h"

using namespace std::string_literals;

namespace BMF {
    struct infer_aux_t {
        FILE *fqh;
        samFile *in;
        samFile *out;
        int cmpkey; // 0 for pos, 1 for unclipped start position
        int mmlim; // Mismatch failure threshold.
        int is_se;
        bam_hdr_t *hdr; // BAM header
        std::unordered_map<std::string, std::string> realign_pairs;
    };
    // Once a set of reads have been selected, their data is put into this struct for collapsing.
    class BamFisherSet {
    const int32_t len;
    std::vector<uint32_t> phred_sums; // Length: 5 * readlen
    std::vector<uint32_t> votes; // Length: 5 * readlen
    std::string name; // Read name
    std::string max_observed_phreds; // Held in memory -33, write out as readable string.
    uint32_t n:31;
    uint32_t is_read1:1;
    public:
        void update_qual(uint8_t *qual) {
            for(int i = 0; i < len; ++i) {
                if(max_observed_phreds[i] < qual[i])
                    max_observed_phreds[i] = qual[i];
            }
        }
        uint32_t get_is_read1() {
            return is_read1;
        }
        BamFisherSet(bam1_t *b) :
        len(b->core.l_qseq),
        phred_sums(len * 5),
        name(bam_get_qname(b)),
        max_observed_phreds(0, len),
        n(1),
        is_read1(!!(b->core.flag & BAM_FREAD1))
        {
            update_qual(bam_get_qual(b));
        }
        void add_rec(bam1_t *b) {
            assert(b->core.l_qseq == len);
            uint8_t *seq = bam_get_seq(b);
            uint8_t *qual = bam_get_qual(b);
            update_qual(qual);
            for(int i = 0; i < len; ++i) {
                const int ind = i * 5 + seqnt2num(bam_seqi(seq, i));
                phred_sums[ind] += qual[i];
                ++votes[ind];
            }
            ++n;
        }
        std::string to_fastq() {
            int i;
            std::stringstream stream;
            stream << '@' << name << " PV:B:I";
            std::string seq(0, len);
            std::vector<uint32_t> agrees(len);
            std::vector<uint32_t> full_quals(len); // igamc calculated
            for(i = 0; i < len; ++i) {
                const int argmaxret = arr_max_u32(phred_sums.data(), i); // 0,1,2,3,4/A,C,G,T,N
                agrees[i] = votes[i * 5 + argmaxret];
                full_quals[i] = pvalue_to_phred(igamc_pvalues(n, LOG10_TO_CHI2(phred_sums[i * 5 + argmaxret])));
                // Mask unconfident base calls
                if(full_quals[i] < 2 || (double)agrees[i] / n < MIN_FRAC_AGREED) {
                    seq[i] = 'N';
                    max_observed_phreds[i] = 2;
                } else seq[i] = num2nuc(argmaxret);
            }
            for(auto pv: full_quals) stream << ',' << pv;
            stream << "\tFA:B:I";
            for(auto agree: agrees) stream << ',' << agree;
            stream << "\tFM:i:" << n << '\n';
            for(auto& phred: max_observed_phreds)
                phred += 33;
            stream << seq << "\n+\n" << max_observed_phreds << '\n';
            return stream.str();
        }
        std::string get_name() {
            return name;
        }
    };
    class BamFisherKing {
    std::unordered_map<int32_t, BamFisherSet> sets;
    public:
        BamFisherKing(dlib::tmp_stack_t *stack) {
            // Name-sort for balanced pairs.
            std::sort(stack->a, stack->a + stack->n, [](const bam1_t *a, const bam1_t *b) {
                return strcmp(bam_get_qname(a), bam_get_qname(b)) < 0;
            });
            std::for_each(stack->a, stack->a + stack->n, [&](bam1_t *b) {
                auto found = sets.find(b->core.l_qseq);
                if(found == sets.end()) {
                    sets.emplace(b->core.l_qseq, BamFisherSet(b));
                } else {
                    found->second.add_rec(b);
                }
            });
        }
        void add_to_hash(infer_aux_t *settings) {
            for(auto set: sets) {
                auto found = settings->realign_pairs.find(set.second.get_name());
                if(found == settings->realign_pairs.end())
                    settings->realign_pairs.emplace(set.second.get_name(), set.second.to_fastq());
                else {
                    if(set.second.get_is_read1()) {
                        fputs(set.second.to_fastq().c_str(), settings->fqh);
                        fputs(found->second.c_str(), settings->fqh);
                    } else {
                        fputs(found->second.c_str(), settings->fqh);
                        fputs(set.second.to_fastq().c_str(), settings->fqh);
                    }
                }
            }
        }
    };
    //extern std::string bam2cppstr(bam1_t *b);
    // In this prototype, we're ignoring the alignment stop, though it should likely be expanded to include it.
} /* namespace BMF */

#endif /* BMF_INFER_H */
