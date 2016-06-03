#include "bmf_main.h"

static int bmftools_usage(int rc)
{
    fprintf(stderr,
                    "Usage: bmftools <subcommand>. See subcommand menus for usage.\n"
                    "-v/--version:            Print bmftools version and exit.\n"
                    "cap:                     Modifies the quality string as function of family metadata.\n"
                    "depth:                   Calculates depth of coverage over a set of bed intervals.\n"
                    "dmp:                     Demultiplex inline barcoded experiments.\n"
                    "err:                     Calculate error rates based on cycle, base call, and quality score.\n"
                    "famstats:                Calculate family size statistics for a bam alignment file.\n"
                    "filter:                  Filter or split a bam file by a set of filters.\n"
                    "inmem:                   Performs dmp fully in memory. RAM-hungry but fast!\n"
                    //"hashdmp:                 Demultiplex inline barcoded experiments that have already been marked.\n"
                    "mark:                    Add tags including unclipped start positions.\n"
                    "rsq:                     Rescue reads with using positional inference to collapse to unique observations in spite of errors in the barcode sequence.\n"
                    "sdmp:                    Demultiplex secondary-index barcoded experiments.\n"
                    "sort:                    Sort for bam rescue.\n"
                    "stack:                   A maximally-permissive yet statistically-thorough variant caller using molecular barcode metadata.\n"
                    "target:                  Calculates on-target rate.\n"
                    "vet:                     Curate variant calls from another variant caller (.bcf) and an indexed alignment file.\n"
            );
    exit(rc);
}

int main(int argc, char *argv[])
{
    if(argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
        return bmftools_usage(EXIT_FAILURE);
    if(strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
        fprintf(stderr,"bmftools version: '%s'.\n", BMF_VERSION);
        exit(EXIT_SUCCESS);
    }
    if(strcmp(argv[1], "sort") == 0) return sort_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "collapse") == 0) return bmf::dmp_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "sdmp") == 0) return bmf::sdmp_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "rsq") == 0) return bmf::rsq_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "hashdmp") == 0) return bmf::hashdmp_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "inmem") == 0) return bmf::hashdmp_inmem_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "famstats") == 0) return bmf::famstats_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "vet") == 0) return bmf::vet_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "err") == 0) return bmf::err_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "mark") == 0) return bmf::mark_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "cap") == 0) return bmf::cap_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "target") == 0) return bmf::target_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "depth") == 0) return bmf::depth_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "stack") == 0) return bmf::stack_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "filter") == 0) return bmf::filter_main(argc - 1, argv + 1);
    fprintf(stderr, "Unrecognized command %s. Abort!\n", argv[1]);
    return EXIT_FAILURE;
}
