#ifndef BMF_VETTER_H
#define BMF_VETTER_H

#include <cstdint>
#include "htslib/khash.h"
#include "lib/stack.h"

namespace bmf {

    int vet_main(int, char **);

    uint64_t NUM_PREALLOCATED_ALLELES = 4uL;
    const char *bmf_header_lines[] =  {
            "##INFO=<ID=BMF_VET,Number=A,Type=Integer,Description=\"1 if the variant passes vetting, 0 otherwise.\">",
            "##INFO=<ID=BMF_QSS,Number=A,Type=Integer,Description=\"Q Score Sum for BMF-passing observations for allele.\">",
            "##INFO=<ID=BMF_UNIOBS,Number=A,Type=Integer,Description=\"Number of unique observations supporting variant at position.\">",
            "##INFO=<ID=BMF_QUANT,Number=A,Type=Integer,Description=\"Estimated quantitation for variant allele.\">",
            "##INFO=<ID=BMF_DUPLEX,Number=A,Type=Integer,Description=\"Number of duplex reads supporting variant at position.\">",
            "##INFO=<ID=BMF_FAIL,Number=A,Type=Integer,Description=\"Number of unique observations at position failing filters.\">",
            "##INFO=<ID=DUPLEX_DEPTH,Number=1,Type=Integer,Description=\"Number of duplex reads passing filters at position.\">",
            "##INFO=<ID=DISC_OVERLAP,Number=1,Type=Integer,Description=\"Number of read pairs at position with discordant base calls.\">",
            "##INFO=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Number of overlapping read pairs combined into single observations at position.\">"
    };

}

#endif /* BMF_VETTER_H */
