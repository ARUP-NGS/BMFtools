#include "bmf_vetter.h"


void vetter_error(char *message, int retcode)
{
	fprintf(stderr, message);
	exit(retcode);
}

void vetter_usage(int retcode)
{
	vetter_error("Not written. Eh.\n", retcode);
}

void vs_open(vetter_settings_t *settings) {
    if((settings->bh = bam_hdr_read(settings->bam)) == NULL)
    	vetter_error("Could not read header from bam. Abort!\n", EXIT_FAILURE);
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

int vetter_main() {
	return 0;
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
	settings->bam = sam_open(settings->bam_path, bam_rmode);

	// Open handles
	vs_open(settings);

	int rc;
	if((rc = vetter_main())) {
		fprintf(stderr, "vetter main failed with retcode %i. Abort mission!\n", rc);
		exit(EXIT_FAILURE);
	}
    // Clean up
    vs_destroy(settings);
    free(settings);
    return 0;
}

int bmf_vetter_main(int argc, char *argv[])
{
	const struct option lopts[] = {VETTER_OPTIONS};
	char bam_rmode[2] = "b";
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
    while ((c = getopt_long(argc, argv, "a:s:m:p:f:b:v:o:?h", lopts, NULL)) >= 0) {
        switch (c) {
        case 'a': params.minFA = strtoul(optarg, NULL, 0); break;
        case 's': params.minFM = strtoul(optarg, NULL, 0); break;
        case 'm': params.minMQ = strtoul(optarg, NULL, 0); break;
        case 'p': params.minPV = strtoul(optarg, NULL, 0); break;
        case 'f': params.minFR = atof(optarg); break;
        case 'b': strcpy(bed, optarg); break;
        case 'v': strcpy(invcf, optarg); break;
        case 'o': strcpy(outvcf, optarg); break;
        case 'h': /* fall-through */
        case '?': vetter_usage(EXIT_SUCCESS);
        default: vetter_error("Unrecognized option. Abort!\n", EXIT_FAILURE);
        }
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
    if(optind >= argc) {
    	vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
    }
    strcpy(inbam, argv[optind]);
    bmf_vetter_core(invcf, inbam, outvcf, bed, bam_rmode, vcf_rmode, vcf_wmode, &params);
    return 0;
}
