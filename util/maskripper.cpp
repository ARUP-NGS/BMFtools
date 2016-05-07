#include <getopt.h>
#include "dlib/bam_util.h"

int usage(char **argv, int retcode=EXIT_FAILURE) {
    fprintf(stderr, "%s <-l output_compression_level> in.bam out.bam\n"
                    "Use - for stdin or stdout.\n", argv[0]);
    return retcode;
}
struct opts_t {
    uint32_t min_trimmed_len;
}

int trim_ns(bam1_t *b, void *data) {
    // Get #Ns at the beginning 
    int i;
    uint8_t *const seq(bam_get_seq(b));
    uint32_t *const cigar(bam_get_cigar(b));

    for(i = 0; bam_seqi(seq, i) == dlib::htseq::HTS_N; ++i);
    const int n_start(i);

    if(i == b->core.l_seq - 1) // all bases are N -- garbage read
         return 1;

    // Get #Ns at the end
    for(i = b->core.l_qseq - 1; bam_seqi(seq, i) == dlib::htseq::HTS_N; --i);
    const int n_end(b->core.l_qseq - 1 - i);

    // Get new length for read
    const int final_len(b->core.l_qseq - n_end - n_start);
    if(final_len < ((opts_t *)data)->min_trimmed_len) // Too short.
        return 1;

    // Get new n_cigar. 

    // If the new n_cigar is changed, up
}

int main(int argc, char *argv[]) {
    if(argc < 3) {
        return usage(argv);
    }
    if(strcmp(argv[1], "--help") == 0) {
        return usage(argv, EXIT_SUCCESS);
    }
    int c;
    char out_mode[4] = "wb";
    opts_t opts{0};
    while((c = getopt(argc, argv, "m:l:h?")) > -1) {
        switch(c) {
        case 'm': min_trimmed_len = (uint32_t)optarg; break;
        case 'l': out_mode[2] = atoi(optarg) % 10 + '0'; break;
        case 'h': case '?': return usage(argv, EXIT_SUCCESS);
        }
    }
    if(argc - 2 != optind)
        LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");

    // Actually this function. You can't really apply a null function....
    void *data = NULL;
    dlib::BamHandle inHandle(argv[optind]);
    dlib::BamHandle outHandle(argv[optind + 1], inHandle.header, "wb");
    // int abstract_single_iter(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux_fn function, void *aux);
    dlib::abstract_single_iter(inHandle.fp, inHandle.header, outHandle.fp,
                               &trim_ns, nullptr);
    return EXIT_SUCCESS;
}
