#include "bmf_vetter.h"


int vetter_usage()
{
	fprintf(stderr, "Not written. Eh.\n");
	exit(0);
}

int bmf_vetter_core(vetter_settings_t *settings, const char *bam_rmode, const char *vcf_rmode,
					const char *vcf_wmode)
{
	vetter_settings_t *settings = (vetter_settings_t *)calloc(1, sizeof(vetter_settings_t));
}

int bmf_vetter_main(int argc, char *argv[])
{
	char bam_rmode[2] = "b";
	char vcf_rmode[4] = "";
	char vcf_wmode[4] = "w";
	char *bam_path = NULL;
	char *vcf_path = NULL;
    while ((c = getopt_long(argc, argv, "o:b:?h", lopts, NULL)) >= 0) {
        switch (c) {
        case 'c': vcf_wmode = "wb"; break;
        case 'b': strcpy(settings->bam_path, optarg); break;
        case 'o': strcpy(settings->out_vcf_path, optarg); break;
        case 'h': /* fall-through */
        case '?': return vetter_usage();
        }
    }
    if(!settings->out_vcf_path) {
    	fprintf(stderr, "[%s] Emitting to stdout.\n");
    	strcpy(settings->out_vcf_path, "-");
    }
    else if(strcmp(strrchr(settings->out_vcf_path, '.'), ".vcf")) == 0 &&
    		vcf_rmode[0] == '\0')
    	strcpy(vcf_rmode, "b");
    else
    	strcpy(vcf_rmode, "rb");
    strcpy(vcf_wmode, strcmp(strrchr(settings->in_vcf_path, '.'), ".vcf") ?
    		"wb": "w");
    if(optind >= argc) {
    	fprintf(stderr, "Insufficient arguments. Input vcf required!\n");
    	vetter_usage();
    }
    strcpy(settings->in_vcf_path, argv[optind]);
    settings->vin = vcf_open(settings->in_vcf_path, vcf_rmode)
    settings->vout = vcf_open(settings->out_vcf, vcf_wmode);
}
