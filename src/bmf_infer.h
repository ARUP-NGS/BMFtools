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

    // Once a set of reads have been selected, their data is put into this class for collapsing.
    class BamFisherSet {
    const int32_t len;
    std::vector<uint32_t> phred_sums; // Length: 5 * readlen
    std::vector<uint32_t> votes; // Length: 5 * readlen
    std::string max_observed_phreds; // Held in memory -33, write out as readable string.
    std::string name;
    uint32_t n:31;
    uint32_t is_read1:1;
    public:
        void update_qual(uint8_t *qual) {
            for(int i = 0; i < len; ++i)
                if(max_observed_phreds[i] < qual[i])
                    max_observed_phreds[i] = qual[i];
        }
        uint32_t get_is_read1() {
            return is_read1;
        }
        std::string const& get_name() {
            return name;
        }
        void make_name(bam1_t *b) {
            if(is_read1) {
                sprintf((char *)name.data(), "collapsed:%i:%i:%i:%i:%i:%i:%i:%i",
                        bam_itag(b, "SU"), bam_itag(b, "MU"), // Unclipped starts for self and mate
                        b->core.tid, b->core.mtid, // Contigs
                        !!(b->core.flag & (BAM_FREVERSE)), !!(b->core.flag & (BAM_FMREVERSE)), // Strandedness combinations
                        bam_itag(b, "LR"), bam_itag(b, "LM") // Read length of self and mate.
                        );
            } else {
                sprintf((char *)name.data(), "collapsed:%i:%i:%i:%i:%i:%i:%i:%i",
                        bam_itag(b, "MU"), bam_itag(b, "SU"),
                        b->core.mtid, b->core.tid,
                        !!(b->core.flag & (BAM_FMREVERSE)), !!(b->core.flag & (BAM_FREVERSE)),
                        bam_itag(b, "LM"), bam_itag(b, "LR")
                        );
            }
            name.resize(strlen(name.data()));
        }
        BamFisherSet(bam1_t *b) :
        len(b->core.l_qseq),
        phred_sums(len * 5),
        max_observed_phreds(0, len),
        name('\0', 60uL),
        n(1),
        is_read1(!!(b->core.flag & BAM_FREAD1))
        {
            update_qual(bam_get_qual(b));
            make_name(b);

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
            stream << '@' << get_name() << " PV:B:I";
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
    };
    class BamFisherKing {
    std::unordered_map<int32_t, BamFisherSet> sets;
    public:
        BamFisherKing(dlib::tmp_stack_t *stack) {
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
                    assert(memcmp(found->second.c_str() + 1,
                                  set.second.get_name().c_str(),
                                  set.second.get_name().size()) == 0);
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
    // In this prototype, we're ignoring the alignment stop, though it should likely be expanded to include it.
} /* namespace BMF */

#endif /* BMF_INFER_H */
