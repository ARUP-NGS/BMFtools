#include "kingfisher.h"

#define dmp_pos(kfp, bufs, argmaxret, i)\
	bufs->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[i * 5 + argmaxret]))));\
	bufs->agrees[i] = kfp->nuc_counts[i * 5 + argmaxret];\
	bufs->cons_seq_buffer[i] = (bufs->cons_quals[i] > 2) ? ARRG_MAX_TO_NUC(argmaxret): 'N';\
    if(bufs->cons_seq_buffer[i] == 'N' && bufs->cons_quals[i] > 2) bufs->cons_quals[i] = 2

void dmp_process_write(KingFisher_t *kfp, FILE *handle, tmpbuffers_t *bufs)
{
	//1. Argmax on the phred_sums arrays, using that to fill in the new seq and
	//buffer[0] = '@'; Set this later?
	bufs->cons_seq_buffer[kfp->readlen] = '\0'; // Null-terminal cons_seq.
	for(int i = 0; i < kfp->readlen; ++i) {
		const int argmaxret = ARRG_MAX(kfp, i);
		dmp_pos(kfp, bufs, argmaxret, i);
	}
		/*
		 * Should I add this back in?
		if(bufs->cons_quals[i] > 1073741824) { // Underflow!
			fprintf(stderr, "Note: phred_sums in underflow: %" PRIu32 ".\n", kfp->phred_sums[i * 4 + argmaxret]);
			bufs->cons_quals[i] = 3114;
		}
		*/
	fill_fa_buffer(kfp, bufs->agrees, bufs->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(kfp, bufs->cons_quals, bufs->PVBuffer);
#if PUTC
	fputc('@', handle);
	fputs(kfp->barcode,handle);
	fputc(' ', handle);
	fputs(bufs->FABuffer, handle);
	fputc('\t', handle);
	fputs(bufs->PVBuffer, handle);
	fputs("\tFP:i:", handle);
	fputc(kfp->pass_fail, handle);
	fputs("\tRV:i:", handle);
	fputc(kfp->n_rc + '0', handle);
	fputs("\tFM:i:", handle);
	fputc(kfp->length + '0', handle);
	fputc('\n', handle);
	fputs(bufs->cons_seq_buffer, handle);
	fputs("\n+\n", handle);
	fputs(kfp->max_phreds, handle);
	fputc('\n', handle);
#else
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tRV:i:%i\tFM:i:%i\n%s\n+\n%s\n", kfp->barcode,
			bufs->FABuffer, bufs->PVBuffer,
			kfp->pass_fail, kfp->n_rc, kfp->length,
			bufs->cons_seq_buffer, kfp->max_phreds);
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
