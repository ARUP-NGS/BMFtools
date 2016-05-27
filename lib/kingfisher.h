#ifndef KINGFISHER_H
#define KINGFISHER_H
#include <assert.h>
#include <math.h>
#include <zlib.h>
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "dlib/cstr_util.h"
#include "include/igamc_cephes.h"
#include "lib/splitter.h"

#ifndef MAX_PV
#    define MAX_PV 3117 // Maximum seen with doubles
#endif
#define HASH_DMP_OFFSET 14
#define FP_OFFSET 9

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif
namespace bmf {

    const double MIN_FRAC_AGREED = 0.5; // Minimum fraction of bases agreed in a family to not "N" the base.

    struct tmpbuffers_t {
        char name_buffer[120];
        char PVBuffer[1000];
        char FABuffer[1000];
        char cons_seq_buffer[SEQBUF_SIZE];
        uint32_t cons_quals[SEQBUF_SIZE];
        uint16_t agrees[SEQBUF_SIZE];
    };


    struct tmpvars_t {
        char *bs_ptr;
        int blen;
        int readlen;
        char key[MAX_BARCODE_LENGTH + 1];
        int l; // For holding ret value for seq.
        tmpbuffers_t *buffers;
    };


    struct kingfisher_t {
        uint16_t *nuc_counts; // Count of nucleotides of this form
        uint32_t *phred_sums; // Sums of -10log10(p-value)
        char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
        int length; // Number of reads in family
        int readlen; // Length of reads
        char barcode[MAX_BARCODE_LENGTH + 1];
        char pass_fail;
    };


    void zstranded_process_write(kingfisher_t *kfpf, kingfisher_t *kfpr, kstring_t *ks, tmpbuffers_t *bufs);
    void dmp_process_write(kingfisher_t *kfp, kstring_t *ks, tmpbuffers_t *bufs, int is_rev);
    int kf_hamming(kingfisher_t *kf1, kingfisher_t *kf2);

    static inline void kfill_both(int readlen, uint16_t *agrees, uint32_t *quals, kstring_t *ks)
    {
        int i;
        kputsnl("FA:B:I", ks);
        for(i = 0; i < readlen; ++i) ksprintf(ks, ",%u", agrees[i]);
        kputsnl("\tPV:B:I", ks);
        for(i = 0; i < readlen; ++i) ksprintf(ks, ",%u", quals[i]);
    }

    static inline void pb_pos(kingfisher_t *kfp, kseq_t *seq, int i) {
        const uint32_t posdata = nuc2num(seq->seq.s[i]) + i * 5;
        ++kfp->nuc_counts[posdata];
        kfp->phred_sums[posdata] += seq->qual.s[i] - 33;
        if(seq->qual.s[i] > kfp->max_phreds[posdata]) kfp->max_phreds[posdata] = seq->qual.s[i];
    }

    static inline void pushback_inmem(kingfisher_t *kfp, kseq_t *seq, int offset, int pass) {
        if(!kfp->length++) {
            kfp->pass_fail = pass + '0';
        } else {
            if(kfp->readlen + offset != (int64_t)seq->seq.l) {
                if(pass) return; // Don't bother, it's an error.
                offset = seq->seq.l - kfp->readlen;
            }
        }
        uint32_t posdata, i;
        for(i = offset; i < seq->seq.l; ++i) {
            posdata = nuc2num(seq->seq.s[i]) + (i - offset) * 5;
            assert(posdata < (unsigned)kfp->readlen * 5);
            ++kfp->nuc_counts[posdata];
            kfp->phred_sums[posdata] += seq->qual.s[i] - 33;
            if(seq->qual.s[i] > kfp->max_phreds[posdata])
                kfp->max_phreds[posdata] = seq->qual.s[i];
        }
    }

    static inline void pushback_kseq(kingfisher_t *kfp, kseq_t *seq, int blen)
    {
        if(!kfp->length++) { // Increment while checking
            kfp->pass_fail = seq->comment.s[FP_OFFSET];
            memcpy(kfp->barcode, seq->comment.s + HASH_DMP_OFFSET, blen);
            kfp->barcode[blen] = '\0';
        }
        for(int i = 0; i < kfp->readlen; ++i) pb_pos(kfp, seq, i);
    }


    /*
     * @func arr_max_u32
     * :param: arr [uint32_t *] 2-d array of values. 5 * index + basecall is the index to use.
     * :param: index [int] Base in read to find the maximum value for.
     * :returns: [int] the nucleotide number for the maximum value at this index in the read.
     */
    CONST static inline int arr_max_u32(uint32_t *arr, int index)
    {
        arr += index * 5;
        return (arr[0] > arr[1]) ? ((arr[0] > arr[2]) ? ((arr[0] > arr[3]) ? (arr[0] > arr[4] ? 0: 4)
                                                                           : (arr[3] > arr[4] ? 3: 4))
                                                      : (arr[2] > arr[3])  ? (arr[2] > arr[4] ? 2: 4)
                                                                           : (arr[3] > arr[4] ? 3: 4))
                                 : ((arr[1] > arr[2]) ? ((arr[1] > arr[3]) ? (arr[1] > arr[4] ? 1: 4)
                                                                           : (arr[3] > arr[4] ? 3: 4))
                                                      : ((arr[2] > arr[3]) ? (arr[2] > arr[4] ? 2: 4)
                                                                           : (arr[3] > arr[4] ? 3: 4)));

    }


    CONST static inline int kfp_argmax(kingfisher_t *kfp, int index)
    {
        return arr_max_u32(kfp->phred_sums, index);
    }

    std::vector<double> get_igamc_threshold(int family_size, int max_phred=MAX_PV, double delta=0.002);
    std::vector<std::vector<double>> get_igamc_thresholds(int max_family_size, int max_phred=MAX_PV, double delta=0.002);

} /* namespace bmf */


#endif /*KINGFISHER_H*/
