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
#include "dlib/cstr_util.h"
#include "dlib/sort_util.h"
#include "dlib/bam_util.h"
#include "bmf_rsq.h"

using namespace std::string_literals;

namespace BMF {
    // Once a set of reads have been selected, their data is put into this struct for collapsing.
    class BamFisherSet {
    std::vector<uint32_t> phred_sums; // Length: 5 * readlen
    std::string name; // Read name
    size_t n;
    public:
        BamFisherSet(bam1_t *b) :
        phred_sums(b->core.l_qseq * 5),
        name(bam_get_qname(b)),
        n(0)
        {
        }
        void add_rec(bam1_t *b) {
            assert((uint32_t)b->core.l_qseq == phred_sums.size() / 5);
            uint8_t *seq = bam_get_seq(b);
            uint8_t *qual = bam_get_qual(b);
            for(int i = 0; i < b->core.l_qseq; ++i)
                phred_sums[i * 5 + seqnt2num(bam_seqi(seq, i))] += qual[i];
            ++n;
        }
        std::string to_fastq() {
            // TODO: WRITE THIS FUNCTION!
            return ""s;
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
        void add_to_hash(std::unordered_map<std::string, std::string>& realign_pairs) {
            for(auto set: sets) {
                realign_pairs.emplace(set.second.get_name(), set.second.to_fastq());
            }
        }
    };
    //extern std::string bam2cppstr(bam1_t *b);
    // In this prototype, we're ignoring the alignment stop, though it should likely be expanded to include it.
} /* namespace BMF */

#endif /* BMF_INFER_H */
