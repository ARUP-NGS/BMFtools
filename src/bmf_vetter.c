#include "bmf_vetter.h"

int max_depth = 20000;


void vetter_error(char *message, int retcode)
{
	fprintf(stderr, message);
	exit(retcode);
}

void vetter_usage(int retcode)
{
	char buf[2000];
	sprintf(buf, "Usage:\nbmftools vet -o <out.vcf [stdout]> <in.vcf> <in.srt.indexed.bam>\n"
			 "Optional arguments:\n"
			 "-b\tPath to bed file to only validate variants in said region\n"
			 "-f\tMinimum fraction of reads in a family agreed on a base call\n"
			 "-a\tMinimum number of reads in a family agreed on a base call\n"
			 "-s\tMinimum number of reads in a family to include a that collapsed observation\n"
			 "-f\tMinimum fraction of reads in a family agreed on that base\n"
			 "-m\tMinimum mapping quality for reads for inclusion\n"
			 "-p\tMinimum calculated p-value on a base call in phred space\n");
	vetter_error(buf, retcode);
}

void vs_open(vetter_settings_t *settings) {
	fprintf(stderr, "Reading header from %s and putting it into a header.\n", settings->bam->fn);
	if((settings->bh = sam_hdr_read(settings->bam)) == NULL)
		vetter_error("Could not read header from bam. Abort!\n", EXIT_FAILURE);
	fprintf(stderr, "Num targets: %.i\n", settings->bh->n_targets);
	settings->vh = bcf_hdr_read(settings->vin);
	bcf_hdr_write(settings->vout, settings->vh);
	if((settings->fai = fai_load(settings->ref_path)) == NULL)
		vetter_error("Could not read Fasta index Abort!\n", EXIT_FAILURE);
	settings->bed = bed_read(settings->bed_path);
}

void vs_destroy(vetter_settings_t *settings) {
	if(hts_close(settings->vin))
		vetter_error("Could not close input vcf. ??? Abort!\n", EXIT_FAILURE);
	if(hts_close(settings->vout))
		vetter_error("Could not close output vcf. ??? Abort!\n", EXIT_FAILURE);
	bam_hdr_destroy(settings->bh);
	bcf_hdr_destroy(settings->vh);
	fai_destroy(settings->fai);
	if(hts_close(settings->bam))
		vetter_error("Could not close input bam. ??? Abort!\n", EXIT_FAILURE);
	bed_destroy(settings->bed);
	free(settings);
	settings = NULL;
}

int bmf_vetter_core(char *invcf, char *inbam, char *outvcf, char *bed,
					const char *bam_rmode, const char *vcf_rmode,
					const char *vcf_wmode, vparams_t *params)
{
	vetter_settings_t *settings = (vetter_settings_t *)calloc(1, sizeof(vetter_settings_t));
	memcpy(&settings->params, params, sizeof(vparams_t));
	strcpy(settings->in_vcf_path, invcf);
	strcpy(settings->bam_path, inbam);
	strcpy(settings->out_vcf_path, outvcf);
	strcpy(settings->bed_path, bed);
	settings->vin = vcf_open(settings->in_vcf_path, vcf_rmode);
	settings->vout = vcf_open(settings->out_vcf_path, vcf_wmode);
#if !NDEBUG
	fprintf(stderr, "bam_rmode: %s. bam_path: %s.\n", bam_rmode, settings->bam_path);
#endif
	htsFormat open_fmt;
	memset(&open_fmt, 0, sizeof(htsFormat));
	open_fmt.category = sequence_data;
	open_fmt.format = bam;
	open_fmt.version.major = 1;
	open_fmt.version.minor = 2;
	settings->bam = sam_open_format(settings->bam_path, "r", &open_fmt);
	if(settings->bam == NULL) {
		fprintf(stderr, "Failed to open input file. Abort mission!");
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "Path to input bam: %s.\n", settings->bam->fn);

	// Open handles
	vs_open(settings);

	// Clean up
	vs_destroy(settings);
	free(settings);
	return 0;
}

int bmf_vetter_main(int argc, char *argv[])
{
	const struct option lopts[] = {VETTER_OPTIONS};
	char bam_rmode[3] = "rb";
	char vcf_rmode[4] = "";
	char vcf_wmode[4] = "w";
	char invcf[200] = "";
	char outvcf[200] = "-";
	char inbam[200] = "";
	char bed[200] = "";
	vparams_t params = {
			.minFA = 0u,
			.minPV = 0u,
			.minFM = 0u,
			.minMQ = 0u,
			.minFR = 0.
	};
	int c;
	while ((c = getopt_long(argc, argv, "d:a:s:m:p:f:b:v:o:?h", lopts, NULL)) >= 0) {
		switch (c) {
		case 'a': params.minFA = strtoul(optarg, NULL, 0); break;
		case 's': params.minFM = strtoul(optarg, NULL, 0); break;
		case 'm': params.minMQ = strtoul(optarg, NULL, 0); break;
		case 'p': params.minPV = strtoul(optarg, NULL, 0); break;
		case 'd': max_depth = atoi(optarg); break;
		case 'f': params.minFR = atof(optarg); break;
		case 'b': strcpy(bed, optarg); break;
		case 'o': strcpy(outvcf, optarg); break;
		case 'h': /* fall-through */
		case '?': vetter_usage(EXIT_SUCCESS);
		default: vetter_error("Unrecognized option. Abort!\n", EXIT_FAILURE);
		}
	}

	if(argc < 3) {
		fprintf(stderr, "Insufficient arguments. Abort!\n");
		vetter_usage(EXIT_FAILURE);
	}
	if(outvcf[0] == '-') {
		fprintf(stderr, "[%s] Emitting to stdout as vcf.\n", __func__);
		strcpy(vcf_wmode, "w");
	}
	else if(strrchr(outvcf, '.') && strcmp(strrchr(outvcf, '.'), ".bcf") == 0 &&
			vcf_rmode[0] == '\0')
		strcpy(vcf_rmode, "rb");
	if(!vcf_rmode[0])
		strcpy(vcf_rmode, "r");
	if(!bed[0])
		vetter_error("Bed file required.\n", EXIT_FAILURE);
	strcpy(vcf_wmode, strrchr(invcf, '.') && strcmp(strrchr(outvcf, '.'), ".bcf") ?
			"w": "wb");
	if(optind + 1 >= argc) {
		vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
	}
	strcpy(invcf, argv[optind]);
	strcpy(inbam, argv[optind + 1]);
	bmf_vetter_core(invcf, inbam, outvcf, bed, bam_rmode, vcf_rmode, vcf_wmode, &params);
	return 0;
}
