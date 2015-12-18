#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

extern int bam_sort(int argc, char *argv[]);
extern int crms_main(int argc, char *argv[]);
extern int fqms_main(int argc, char *argv[]);
extern int bam_rsq(int argc, char *argv[]);
extern int hash_dmp_main(int argc, char *argv[]);
extern int famstats_main(int argc, char *argv[]);
extern int bmf_vetter_main(int argc, char *argv[]);
extern int err_main(int argc, char *argv[]);
extern int mark_unclipped_main(int argc, char *argv[]);
extern int cap_qscore_main(int argc, char *argv[]);

static int bmftools_usage(int rc)
{
	fprintf(stderr, "Usage: bmftools <subcommand>. See subcommand menus for usage.\n");
	fprintf(stderr, "sort:            Sort for bam rescue.\n"
                    "dmp:             Demultiplex inline barcoded experiments.\n"
                    "sdmp:            Demultiplex secondary-index barcoded experiments.\n"
                    "hashdmp:         Demultiplex inline barcoded experiments that have already been marked.\n"
                    "rsq:             Rescue bmf-sorted or ucs-sorted bam alignments.\n"
                    "err:             Calculate error rates based on cycle, base call, and quality score.\n"
                    "famstats:        Calculate family size statistics for a bam alignment file.\n"
                    "vet:             Curate variant calls from another variant caller (.vcf) and a bam alignment.\n"
                    "mark_unclipped:  Add unclipped start position as annotation for both read and mate.\n"
                    "cap:             Modifies the quality string as function of family metadata.\n"
			);
	exit(rc);
}

int main(int argc, char *argv[])
{
	if(argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
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
		return hash_dmp_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "famstats") == 0)
		return famstats_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "vet") == 0)
		return bmf_vetter_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "err") == 0)
		return err_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "mark_unclipped") == 0)
		return mark_unclipped_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "cap") == 0)
		return cap_qscore_main(argc - 1, argv + 1);
	fprintf(stderr, "Unrecognized command %s. Abort!\n", argv[1]);
	return EXIT_FAILURE;
}
