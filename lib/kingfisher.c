#include "kingfisher.h"

#define dmp_pos(kfp, bufs, argmaxret, i)\
	bufs->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[i * 5 + argmaxret]))));\
	bufs->agrees[i] = kfp->nuc_counts[i * 5 + argmaxret];\
	bufs->cons_seq_buffer[i] = (bufs->cons_quals[i] > 2) ? ARRG_MAX_TO_NUC(argmaxret): 'N';\
    if(bufs->cons_seq_buffer[i] == 'N' && bufs->cons_quals[i] > 2) bufs->cons_quals[i] = 2

// TODO: rewrite dmp_process_write and stranded_process_write to write the PV/FA strings straight to output
// rather than writing to a temporary object and writing that later.

void dmp_process_write(KingFisher_t *kfp, FILE *handle, tmpbuffers_t *bufs, int is_rev)
{
	for(int i = 0; i < kfp->readlen; ++i) {
		const int argmaxret = ARRG_MAX(kfp, i);
		dmp_pos(kfp, bufs, argmaxret, i);
	}
	fill_fa_buffer(kfp->readlen, bufs->agrees, bufs->FABuffer);
	fill_pv_buffer(kfp->readlen, bufs->cons_quals, bufs->PVBuffer);
#if PUTC
	fputc('@', handle);
	fputs(kfp->barcode,handle);
	fputc(' ', handle);
	fputs(bufs->FABuffer, handle);
	fputc('\t', handle);
	fputs(bufs->PVBuffer, handle);
	fputs("\tFP:i:", handle);
	fputc(kfp->pass_fail, handle);
	fprintf(handle, "\tFM:i:%i", kfp->length);
	fprintf(handle, "\tRV:i:%i\n", is_rev ? kfp->length: 0);
	fputs(bufs->cons_seq_buffer, handle);
	fputs("\n+\n", handle);
	fputs(kfp->max_phreds, handle);
	fputc('\n', handle);
#else
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tFM:i:%i\n%s\n+\n%s\n", kfp->barcode,
			bufs->FABuffer, bufs->PVBuffer,
			kfp->pass_fail, kfp->length,
			bufs->cons_seq_buffer, kfp->max_phreds);
#endif
	return;
}

// kfp forward, kfp reverse
void stranded_process_write(KingFisher_t *kfpf, KingFisher_t *kfpr, FILE *handle, tmpbuffers_t *bufs)
{
	for(int i = 0; i < kfpf->readlen; ++i) {
		const int argmaxretf = ARRG_MAX(kfpf, i), argmaxretr = ARRG_MAX(kfpr, i);
        if(argmaxretf == argmaxretr) { // Both strands supported the same base call.
            kfpf->phred_sums[i * 5 + argmaxretf] += kfpr->phred_sums[i * 5 + argmaxretf];
            kfpf->nuc_counts[i * 5 + argmaxretf] += kfpr->nuc_counts[i * 5 + argmaxretf];
            dmp_pos(kfpf, bufs, argmaxretf, i);
            if(kfpr->max_phreds[i] > kfpf->max_phreds[i]) kfpf->max_phreds[i] = kfpr->max_phreds[i];
        }
        else if(argmaxretf == 5) { // Forward is Nd and reverse is not. Reverse call is probably right.
            kfpf->phred_sums[i * 5 + argmaxretr] += kfpr->phred_sums[i * 5 + argmaxretr];
            kfpf->nuc_counts[i * 5 + argmaxretr] += kfpr->nuc_counts[i * 5 + argmaxretr];
            dmp_pos(kfpf, bufs, argmaxretr, i);
            kfpf->max_phreds[i] = kfpr->max_phreds[i];
        }
        else if(argmaxretr == 5) { // Forward is Nd and reverse is not. Reverse call is probably right.
            kfpf->phred_sums[i * 5 + argmaxretf] += kfpr->phred_sums[i * 5 + argmaxretf];
            kfpf->nuc_counts[i * 5 + argmaxretf] += kfpr->nuc_counts[i * 5 + argmaxretf];
            dmp_pos(kfpf, bufs, argmaxretf, i);
            // Don't update max_phreds, since the max phred is already here.
        }
        else
            bufs->cons_quals[i] = 0, bufs->agrees[i] = 0, bufs->cons_seq_buffer[i] = 'N';
	}
	fill_fa_buffer(kfpf->readlen, bufs->agrees, bufs->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(kfpf->readlen, bufs->cons_quals, bufs->PVBuffer);
#if PUTC
	fputc('@', handle);
	fputs(kfpf->barcode, handle);
	fputc(' ', handle);
	fputs(bufs->FABuffer, handle);
	fputc('\t', handle);
	fputs(bufs->PVBuffer, handle);
	fputs("\tFP:i:", handle);
	fputc(kfpf->pass_fail, handle);
    fprintf(handle, "\tRV:i:%i\tFM:i:%i\n", kfpr->length, kfpf->length + kfpr->length);
	fputs(bufs->cons_seq_buffer, handle);
	fputs("\n+\n", handle);
	fputs(kfpf->max_phreds, handle);
	fputc('\n', handle);
#else
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tRV:i:%i\tFM:i:%i\n%s\n+\n%s\n", kfpf->barcode + 1,
			bufs->FABuffer, bufs->PVBuffer,
			kfpf->pass_fail, kfpr->length, kfpf->length + kfpr->length,
			bufs->cons_seq_buffer, kfpf->max_phreds);
#endif
	return;
}

KingFisher_t *init_kfp(size_t readlen)
{
	KingFisher_t *ret = (KingFisher_t *)calloc(1, sizeof(KingFisher_t));
	ret->readlen = readlen;
	ret->max_phreds = (char *)calloc(readlen + 1, sizeof(char)); // Keep track of the maximum phred score observed at position.
	ret->nuc_counts = (uint16_t *)calloc(readlen * 5, sizeof(uint16_t));
	ret->phred_sums = (uint32_t *)calloc(readlen * 5, sizeof(uint32_t));
	ret->pass_fail = '1';
	return ret;
}

void destroy_kf(KingFisher_t *kfp)
{
	cond_free(kfp->nuc_counts);
	cond_free(kfp->phred_sums);
	cond_free(kfp->max_phreds);
	free(kfp);
	kfp = NULL;
}
