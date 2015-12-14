#include "kingfisher.h"

#define dmp_pos(kfp, bufs, argmaxret, i)\
	bufs->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[i * 5 + argmaxret]))));\
	bufs->agrees[i] = kfp->nuc_counts[i * 5 + argmaxret];\
	bufs->cons_seq_buffer[i] = (bufs->cons_quals[i] > 2 && (double)bufs->agrees[i] / kfp->length > MIN_FRAC_AGREED) ? ARRG_MAX_TO_NUC(argmaxret): 'N';\
    if(bufs->cons_seq_buffer[i] == 'N' && bufs->cons_quals[i] > 2) bufs->cons_quals[i] = 2

// TODO: rewrite dmp_process_write and stranded_process_write to write the PV/FA strings straight to output
// rather than writing to a temporary object and writing that later.

void dmp_process_write(KingFisher_t *kfp, FILE *handle, tmpbuffers_t *bufs, int is_rev)
{
	int i;
	for(i = 0; i < kfp->readlen; ++i) {
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
	for(i = 0; i < kfp->length; ++i)
		fputc(kfp->max_phreds[nuc2num(bufs->cons_seq_buffer[i])][i], handle);
	fputc('\n', handle);
#else
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tFM:i:%i\tRV:i:%i\n%s\n+\n", kfp->barcode + 1,
			bufs->FABuffer, bufs->PVBuffer,
			kfp->pass_fail, kfp->length, is_rev ? kfp->length: 0,
			bufs->cons_seq_buffer);
	for(i = 0; i < kfp->length; ++i)
		fputc(kfp->max_phreds[nuc2num(bufs->cons_seq_buffer[i])][i], handle);
	fputc('\n', handle);
#endif
	return;
}

// kfp forward, kfp reverse
void stranded_process_write(KingFisher_t *kfpf, KingFisher_t *kfpr, FILE *handle, tmpbuffers_t *bufs)
{
	int i;
	for(i = 0; i < kfpf->readlen; ++i) {
		const int argmaxretf = ARRG_MAX(kfpf, i), argmaxretr = ARRG_MAX(kfpr, i);
        if(argmaxretf == argmaxretr) { // Both strands supported the same base call.
            kfpf->phred_sums[i * 5 + argmaxretf] += kfpr->phred_sums[i * 5 + argmaxretf];
            kfpf->nuc_counts[i * 5 + argmaxretf] += kfpr->nuc_counts[i * 5 + argmaxretf];
            dmp_pos(kfpf, bufs, argmaxretf, i);
            if(kfpr->max_phreds[argmaxretf][i] > kfpf->max_phreds[argmaxretf][i])
            	kfpf->max_phreds[argmaxretf][i] = kfpr->max_phreds[argmaxretf][i];
        }
        else if(argmaxretf == 4) { // Forward is Nd and reverse is not. Reverse call is probably right.
            kfpf->phred_sums[i * 5 + argmaxretr] += kfpr->phred_sums[i * 5 + argmaxretr];
            kfpf->nuc_counts[i * 5 + argmaxretr] += kfpr->nuc_counts[i * 5 + argmaxretr];
            dmp_pos(kfpf, bufs, argmaxretr, i);
            kfpf->max_phreds[argmaxretr][i] = kfpr->max_phreds[argmaxretr][i];
        }
        else if(argmaxretr == 4) { // Forward is Nd and reverse is not. Reverse call is probably right.
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
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tFM:i:%i\tRV:i:%i\n%s\n+\n", kfpf->barcode + 1,
			bufs->FABuffer, bufs->PVBuffer,
			kfpf->pass_fail, kfpf->length, kfpr->length,
			bufs->cons_seq_buffer);
#if !NDEBUG
	fprintf(stderr, "[D:%s] Quality string I'm about to write out for consensus sequence '%s':\t",
			__func__, bufs->cons_seq_buffer);
	for(i = 0; i < kfpf->length; ++i)
		fputc(kfpf->max_phreds[nuc2num(bufs->cons_seq_buffer[i])][i], stderr);
	fputc('\n', stderr);
	for(i = 0; i < kfpf->length; ++i) {
		fprintf(stderr, "At position %i with base call %c, max_phreds for A is %c.\n", i + 1, bufs->cons_seq_buffer[i],
				kfpf->max_phreds[0][i]);
		fprintf(stderr, "At position %i with base call %c, max_phreds for C is %c.\n", i + 1, bufs->cons_seq_buffer[i],
				kfpf->max_phreds[1][i]);
		fprintf(stderr, "At position %i with base call %c, max_phreds for G is %c.\n", i + 1, bufs->cons_seq_buffer[i],
				kfpf->max_phreds[2][i]);
		fprintf(stderr, "At position %i with base call %c, max_phreds for T is %c.\n", i + 1, bufs->cons_seq_buffer[i],
				kfpf->max_phreds[3][i]);
		fprintf(stderr, "At position %i with base call %c, max_phreds for N is %c.\n", i + 1, bufs->cons_seq_buffer[i],
				kfpf->max_phreds[4][i]);
	}
#endif
	for(i = 0; i < kfpf->length; ++i)
		fputc(kfpf->max_phreds[nuc2num(bufs->cons_seq_buffer[i])][i], handle);
	fputc('\n', handle);
	return;
}

KingFisher_t *init_kfp(size_t readlen)
{
	KingFisher_t *ret = (KingFisher_t *)calloc(1, sizeof(KingFisher_t));
	ret->readlen = readlen;
	ret->max_phreds = (char **)malloc(5 * sizeof(char *));
	for(int i = 0; i < 5; ++i) {
		ret->max_phreds[i] = (char *)calloc(readlen + 1, sizeof(char));
	}
	ret->nuc_counts = (uint16_t *)calloc(readlen * 5, sizeof(uint16_t));
	ret->phred_sums = (uint32_t *)calloc(readlen * 5, sizeof(uint32_t));
	ret->pass_fail = '1';
	return ret;
}

void destroy_kf(KingFisher_t *kfp)
{
	cond_free(kfp->nuc_counts);
	cond_free(kfp->phred_sums);
	for(int i = 0; i < 5; ++i) {
		cond_free(kfp->max_phreds[i]);
	}
	cond_free(kfp->max_phreds);
	free(kfp);
	kfp = NULL;
}
