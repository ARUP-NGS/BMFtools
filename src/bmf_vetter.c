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

int bmf_vetter_core(char *invcf, char *inbam, char *outvcf,
					const char *bam_rmode, const char *vcf_rmode,
					const char *vcf_wmode)
{
	vetter_settings_t *settings = (vetter_settings_t *)calloc(1, sizeof(vetter_settings_t));
	strcpy(settings->in_vcf_path, invcf);
	strcpy(settings->bam_path, inbam);
	strcpy(settings->out_vcf_path, outvcf);

	// Open handles
    settings->vin = vcf_open(settings->in_vcf_path, vcf_rmode);
    settings->vout = vcf_open(settings->out_vcf, vcf_wmode);
    settings->bam = sam_open(settings->bam_path, bam_rmode);
    if((settings->bh = bam_hdr_read(settings->bam)) == NULL)
    	vetter_error("Could not read header from bam. Abort!\n", EXIT_FAILURE);
    settings->vh = bcf_hdr_read(settings->in_vcf_path);
    bcf_hdr_write(settings->vout, settings->vh);
    settings->fai = fai_load(settings->ref);

    // Clean up
    if(hts_close(settings->vin))
    	vetter_error("Could not close input vcf. ??? Abort!\n", EXIT_FAILURE);
    if(hts_close(settings->vout))
    	vetter_error("Could not close output vcf. ??? Abort!\n", EXIT_FAILURE);
    bam_hdr_destroy(settings->bh);
    bcf_hdr_destroy(settings->vh);
    fai_destroy(settings->fai);
    if(hts_close(settings->bam))
    	vetter_error("Could not close input bam. ??? Abort!\n", EXIT_FAILURE);
    free(settings);
}

int bmf_vetter_main(int argc, char *argv[])
{
	const struct option lopts[] = {VETTER_OPTIONS};
	char bam_rmode[2] = "b";
	char vcf_rmode[4] = "";
	char vcf_wmode[4] = "w";
	char invcf[200] = "";
	char outvcf[200] = "-";
	char inbam[200] "";
    while ((c = getopt_long(argc, argv, "v:o:?h", lopts, NULL)) >= 0) {
        switch (c) {
        case 'v': strcpy(invcf, optarg); break;
        case 'o': strcpy(outvcf, optarg); break;
        case 'h': /* fall-through */
        case '?': return vetter_usage(EXIT_SUCCESS);
        default: return vetter_error("Unrecognized option. Abort!\n", EXIT_FAILURE);
        }
    }
    if(outvcf[0] == '-') {
    	fprintf(stderr, "[%s] Emitting to stdout as vcf.\n", __func__);
    	strcpy(vcf_wmode, "w");
    }
    else if(strrchr(outvcf, '.') && strcmp(strrchr(outvcf, '.'), ".bcf")) == 0 &&
    		vcf_rmode[0] == '\0')
    	strcpy(vcf_rmode, "rb");
    if(!vcf_rmode[0])
    	strcpy(vcf_rmode, "r");
    strcpy(vcf_wmode, strrchr(settings->in_vcf_path, '.') && strcmp(strrchr(outvcf), '.'), ".bcf") ?
    		"w": "wb");
    if(optind >= argc) {
    	vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
    }
    strcpy(inbam, argv[optind]);
    bmf_vetter_core(invcf, inbam, outvcf, bam_rmode, vcf_rmode, vcf_wmode);
}
