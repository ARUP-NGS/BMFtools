#include "kingfisher.h"

namespace BMF {

    const double MIN_FRAC_AGREED = 0.5; // Minimum fraction of bases agreed in a family to not "N" the base.

    #define dmp_pos(kfp, bufs, argmaxret, i, index, diffcount)\
        do {\
            bufs->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[index]))));\
            bufs->agrees[i] = kfp->nuc_counts[index];\
            diffcount -= bufs->agrees[i];\
            if(argmaxret != 4) diffcount -= kfp->nuc_counts[i * 5 + 4];\
            if(bufs->cons_quals[i] <= 2) {\
                bufs->cons_quals[i] = 2;\
                bufs->cons_seq_buffer[i] = 'N';\
            } else if((double)bufs->agrees[i] / kfp->length > MIN_FRAC_AGREED) {\
                bufs->cons_seq_buffer[i] = num2nuc(argmaxret);\
            } else {\
                bufs->cons_quals[i] = 2;\
                bufs->cons_seq_buffer[i] = 'N';\
            }\
        } while(0)

    void dmp_process_write(kingfisher_t *kfp, kstring_t *ks, tmpbuffers_t *bufs, int is_rev)
    {
        int i, diffs = kfp->length * kfp->readlen;
        for(i = 0; i < kfp->readlen; ++i) {
            const int argmaxret = kfp_argmax(kfp, i);
            const int index = argmaxret + i * 5;
            dmp_pos(kfp, bufs, argmaxret, i, index, diffs);
        }
        ksprintf(ks, "@%s ", kfp->barcode + 1);
        kfill_both(kfp->readlen, bufs->agrees, bufs->cons_quals, ks);
        ksprintf(ks, "\tFP:i:%c\tFM:i:%i\tRV:i:%i\tNF:f:%0.6f\tDR:i:0\n%s\n+\n",
                kfp->pass_fail, kfp->length, is_rev ? kfp->length: 0,
                (double) diffs / kfp->length, bufs->cons_seq_buffer);
        for(i = 0; i < kfp->readlen; ++i) kputc(kfp->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + 5 * i], ks);
        kputc('\n', ks);
    }

    // kfp forward, kfp reverse
    // Note: You print kfpf->barcode + 1 because that skips the F/R/Z char.
    void zstranded_process_write(kingfisher_t *kfpf, kingfisher_t *kfpr, kstring_t *ks, tmpbuffers_t *bufs)
    {
        const int FM = kfpf->length + kfpr->length;
        int i, diffs = FM * kfpf->readlen;
        int index;
        for(i = 0; i < kfpf->readlen; ++i) {
            const int argmaxretf = kfp_argmax(kfpf, i), argmaxretr = kfp_argmax(kfpr, i);
            if(argmaxretf == argmaxretr) { // Both strands supported the same base call.
                index = i * 5 + argmaxretf;
                kfpf->phred_sums[index] += kfpr->phred_sums[index];
                kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
                dmp_pos(kfpf, bufs, argmaxretf, i, index, diffs);
                if(kfpr->max_phreds[index] > kfpf->max_phreds[index]) kfpf->max_phreds[index] = kfpr->max_phreds[index];
            } else if(argmaxretf == 4) { // Forward is Nd and reverse is not. Reverse call is probably right.
                index = i * 5 + argmaxretr;
                kfpf->phred_sums[index] += kfpr->phred_sums[index];
                kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
                dmp_pos(kfpf, bufs, argmaxretr, i, index, diffs);
                kfpf->max_phreds[index] = kfpr->max_phreds[index];
            } else if(argmaxretr == 4) { // Forward is Nd and reverse is not. Reverse call is probably right.
                index = i * 5 + argmaxretf;
                kfpf->phred_sums[index] += kfpr->phred_sums[index];
                kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
                dmp_pos(kfpf, bufs, argmaxretf, i, index, diffs);
                // Don't update max_phreds, since the max phred is already here.
            } else bufs->cons_quals[i] = 0, bufs->agrees[i] = 0, bufs->cons_seq_buffer[i] = 'N';
        }
        ksprintf(ks, "@%s ", kfpf->barcode + 1);
        // Add read name
        kfill_both(kfpf->readlen, bufs->agrees, bufs->cons_quals, ks);
        ksprintf(ks, "\tFP:i:%c\tFM:i:%i\tRV:i:%i\tNF:f:%f\tDR:i:%i\n%s\n+\n", kfpf->pass_fail,
                FM, kfpr->length, (double) diffs / FM, kfpf->length && kfpr->length,
                bufs->cons_seq_buffer);
        for(i = 0; i < kfpf->readlen; ++i)
            kputc(kfpf->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + 5 * i], ks);
        kputc('\n', ks);
        //const int ND = get_num_differ
        return;
    }

} /* namespace BMF */
