#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

extern int bam_sort(int argc, char *argv[]);
extern int crms_main(int argc, char *argv[]);
extern int fqms_main(int argc, char *argv[]);
extern int bam_rsq(int argc, char *argv[]);
extern int khash_dmp_main(int argc, char *argv[]);
extern int famstats_main(int argc, char *argv[]);
extern int bmf_vetter_main(int argc, char *argv[]);

int bmftools_usage(int rc)
{
	fprintf(stderr, "Usage: bmftools <subcommand>. See subcommand menus for usage.\n");
	fprintf(stderr, "sort:\tSort for bam rescue.\n"
                    "dmp:\tDemultiplex inline barcoded experiments.\n"
                    "sdmp:\tDemultiplex secondary-index barcoded experiments.\n"
                    "rsq:\tRescue bmf-sorted or ucs-sorted bam alignments.\n"
                    "famstats:\tCalculate family size statistics for a bam alignment file.\n"
                    "vet:\tCurate variant calls from another variant caller (.vcf) and a bam alignment.\n"
			);
	exit(rc);
}

int main(int argc, char *argv[])
{
	if(argc == 1)
		return bmftools_usage(EXIT_FAILURE);
	if(strcmp(argv[1], "sort") == 0)
		return bam_sort(argc - 1, argv + 1);
	else if(strcmp(argv[1], "dmp") == 0)
		return crms_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "sdmp") == 0)
		return fqms_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "rsq") == 0)
		return bam_rsq(argc - 1, argv + 1);
	else if(strcmp(argv[1], "hashdmp") == 0)
		return khash_dmp_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "famstats") == 0)
		return famstats_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "vet") == 0)
		return bmf_vetter_main(argc - 1, argv + 1);
	fprintf(stderr, "Unrecognized command %s. Abort!\n", argv[1]);
}
