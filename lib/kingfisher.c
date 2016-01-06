#include "kingfisher.h"

#define dmp_pos(kfp, bufs, argmaxret, i, index, diffcount)\
	do {\
		bufs->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[index]))));\
		bufs->agrees[i] = kfp->nuc_counts[index];\
		diffcount -= bufs->agrees[i];\
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

// TODO: rewrite dmp_process_write and stranded_process_write to write the PV/FA strings straight to output
// rather than writing to a temporary object and writing that later.

void dmp_process_write(KingFisher_t *kfp, FILE *handle, tmpbuffers_t *bufs, int is_rev)
{
	int i, diffs = kfp->length * kfp->readlen;
	for(i = 0; i < kfp->readlen; ++i) {
		const int argmaxret = kfp_argmax(kfp, i);
		const int index = argmaxret + i * 5;
		dmp_pos(kfp, bufs, argmaxret, i, index, diffs);
	}
	fill_fa(kfp->readlen, bufs->agrees, bufs->FABuffer);
	fill_pv(kfp->readlen, bufs->cons_quals, bufs->PVBuffer);
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
	for(i = 0; i < kfp->readlen; ++i)
		fputc(kfp->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + i * 5], handle);
	fputc('\n', handle);
#else
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tFM:i:%i\tRV:i:%i\tNF:f:%0.6f\n%s\n+\n", kfp->barcode + 1,
			bufs->FABuffer, bufs->PVBuffer,
			kfp->pass_fail, kfp->length, is_rev ? kfp->length: 0,
			(double) diffs / kfp->length,
			bufs->cons_seq_buffer);
	for(i = 0; i < kfp->readlen; ++i) fputc(kfp->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + 5 * i], handle);
	fputc('\n', handle);
#endif
	return;
}

// kfp forward, kfp reverse
// Note: You print kfpf->barcode + 1 because that skips the F/R/Z char.
void stranded_process_write(KingFisher_t *kfpf, KingFisher_t *kfpr, FILE *handle, tmpbuffers_t *bufs)
{
	int i, diffs=kfpf->length * kfpf->readlen;
	for(i = 0; i < kfpf->readlen; ++i) {
		const int argmaxretf = kfp_argmax(kfpf, i), argmaxretr = kfp_argmax(kfpr, i);
		if(argmaxretf == argmaxretr) { // Both strands supported the same base call.
			const int index = i * 5 + argmaxretf;
			kfpf->phred_sums[index] += kfpr->phred_sums[index];
			kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
			dmp_pos(kfpf, bufs, argmaxretf, i, index, diffs);
			if(kfpr->max_phreds[index] > kfpf->max_phreds[index]) kfpf->max_phreds[index] = kfpr->max_phreds[index];
		} else if(argmaxretf == 4) { // Forward is Nd and reverse is not. Reverse call is probably right.
			const int index = i * 5 + argmaxretr;
			kfpf->phred_sums[index] += kfpr->phred_sums[index];
			kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
			dmp_pos(kfpf, bufs, argmaxretr, i, index, diffs);
			kfpf->max_phreds[index] = kfpr->max_phreds[index];
		} else if(argmaxretr == 4) { // Forward is Nd and reverse is not. Reverse call is probably right.
			const int index = i * 5 + argmaxretf;
			kfpf->phred_sums[index] += kfpr->phred_sums[index];
			kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
			dmp_pos(kfpf, bufs, argmaxretf, i, index, diffs);
			// Don't update max_phreds, since the max phred is already here.
		} else bufs->cons_quals[i] = 0, bufs->agrees[i] = 0, bufs->cons_seq_buffer[i] = 'N';
	}
	fill_fa(kfpf->readlen, bufs->agrees, bufs->FABuffer);
	fill_pv(kfpf->readlen, bufs->cons_quals, bufs->PVBuffer);
	//const int ND = get_num_differ
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tFM:i:%i\tRV:i:%i\tNF:f:%f\n%s\n+\n", kfpf->barcode + 1,
			bufs->FABuffer, bufs->PVBuffer,
			kfpf->pass_fail, kfpf->length, kfpr->length,
			(double) diffs / kfpf->length,
			bufs->cons_seq_buffer);
	for(i = 0; i < kfpf->readlen; ++i)
		fputc(kfpf->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + 5 * i], handle);
	fputc('\n', handle);
	return;
}
